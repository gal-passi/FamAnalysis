import os
from os.path import basename
import numpy as np
from Bio import SeqIO
import json
import warnings
import glob
import pandas as pd
import pickle
from functools import lru_cache
from pathlib import Path
from utils import afm_range_read, name_for_esm, sequence_from_esm_df, name_for_cpt, \
    sequence_from_cpt_df, seacrh_by_ref_seq, warn_if, esm_seq_logits, esm_process_long_sequences, \
    sequence_from_afm_df
import Connections
from definitions import *

ROOT = os.path.join(os.path.dirname(__file__))


def convert_csv_to_parquet(csv_path, model):
    """Converts CSV to Parquet if it doesn't exist and returns the Parquet file path."""
    filename = os.path.basename(csv_path).replace(".csv", ".parquet")
    if model == 'ESM':
        parquet_path = os.path.join(ESM_PARQUET_DIR, filename)
    elif model == 'EVE':
        parquet_path = os.path.join(EVE_PARQUET_DIR, filename)
    else:
        raise ValueError("Model must be 'ESM' or 'EVE'")

    if not os.path.exists(parquet_path):
        if model == 'EVE':
            df = pd.read_csv(csv_path, low_memory=False)
        elif model == 'ESM':
            df = pd.read_csv(csv_path)
        else:
            raise ValueError("Model must be 'ESM' or 'EVE'")
        df.to_parquet(parquet_path)
        print(f"Converted: {csv_path} → {parquet_path}")

    return parquet_path


@lru_cache(maxsize=200)
def load_parquet_cached(csv_path, model):
    """Loads a Parquet file from disk with LRU caching."""
    parquet_path = convert_csv_to_parquet(csv_path, model)
    if model == 'EVE':
        return pd.read_parquet(parquet_path).fillna(-1)
    elif model == 'ESM':
        return pd.read_parquet(parquet_path)
    else:
        raise ValueError("Model must be 'ESM' or 'EVE'")


class ProteinAnalyzer:
    STRINGS_REVERSE_INDEX_ROOT = f"{ROOT}/DB/strings/strings_rindex.txt"
    PDBS_DB_ROOT = f"{ROOT}/DB/pdbs"

    def __init__(self, *proteins, verbose_level=1):
        self.proteins = {protein.name: protein for protein in proteins}
        self._cpt_ingene = {basename(p)[:-CPT_TRUNCATE_TAIL] for p in glob.glob(pjoin(CPT_INGENE_PATH, '*.csv.gz'))}
        self._cpt_exgene = {basename(p)[:-CPT_TRUNCATE_TAIL] for p in glob.glob(pjoin(CPT_EXGENE_PATH, '*.csv.gz'))}
        self.eve_records = {basename(p)[:-EVE_TRUNCATE_TAIL] for p in glob.glob(pjoin(EVE_VARIANTS_PATH, '*.csv'))}
        self._unip = Connections.Uniport(verbose_level=verbose_level)
        with open(EVE_INDEX_PATH_2, "rb") as fp:
            self.eve_index = pickle.load(fp)
        with open(EVE_INVERSE_INDEX, "rb") as fp:
            self.eve_reverse_index = pickle.load(fp)
        with open(AFM_DIRECTORY_PATH, 'r') as file:
            self.af_index = json.load(file)
        with open(AFM_RANGES_PATH, 'r') as file:
            self.af_ranges = json.load(file)
        with open(ESM_INDEX_PATH, 'r') as file:
            self.esm_index = json.load(file)
        self.pioneer_sequences, self.pioneer_predictions = None, None  # high memory usage load on demand

    def _load_pioneer_data(self):
        with open(PIONEER_SEQ_PATH, 'r') as file:
            self.pioneer_sequences = json.load(file)
        with open(PIONEER_PRED_PATH, 'r') as file:
            self.pioneer_predictions = json.load(file)

    def erase_pioneer_data(self):
        del self.pioneer_sequences
        del self.pioneer_predictions
        self.pioneer_sequences, self.pioneer_predictions = None, None

    @staticmethod
    def _find_chains(chains_dict, mutation):
        """
        return the ids of all chains containing the ref_seq
        :param chains_dict: {chain_id: seq}
        :param mut: Mutation object
        :return: {(chain_ids, location)}
        """
        chains = set()
        for ref in mutation.ref_seqs.values():
            offset = ref.find(mutation.origAA)
            res = {(c_id, seq.find(ref) + offset) for c_id, seq in chains_dict.items() if ref in seq}
            chains = chains.union(res)
        return chains

    @staticmethod
    def pdb_chains(pdb_path):
        """
        :param pdb: path to bdn file
        :param ref_seq: reference sequence of amino acids
        :return: {chain id: Bio.seq}
        """
        res = {record.id: record.seq for record in SeqIO.parse(pdb_path, 'pdb-seqres')}
        return {id.split(':')[1] if ':' in id else id: seq for id, seq in res.items()}

    @staticmethod
    def search_eve_record(prot_name):
        """
        searches for an existing evemodel data
        :param prot_name as in Protein.name:
        :return: bool
        """
        find_file = lambda x: os.path.basename(x)[:-10]
        return prot_name in list(map(find_file, glob.glob(pjoin(EVE_DATA_PATH, '*'))))

    # TODO change eve db to work with unip - search by name might miss some proteins
    def score_mutation_evemodel(self, protein, mutation, prot_name=None, download_records=False, optimized=0):
        """
        looks for an evemodel score for the mutation given
        :param download_records: bool whether to download new eve records. This should be True is user choose setup
               with partial eve data.
        :param protein: Protein object
        :param mutation: Mutation object
        :return: (score, prediction) if not found (-1,-1)
        """
        search_name = protein.name if not prot_name else prot_name
        #  Case 1: Protein reference is directly found in eve records
        if search_name in self.eve_records:
            res = self._eve_interperter(search_name, mutation, optimized)
            if res != (-1, -1):
                return res
        # Case 2: Protein has an alias directly found in eve records
        for alias in protein.aliases:
            if alias in self.eve_records:
                res = self._eve_interperter(alias, mutation, optimized)
                if res != (-1, -1):
                    return res
        #  Case 3: Protein main reference is an alias of a record directly covered by eve
        if search_name in self.eve_reverse_index:
            for name in set(self.eve_reverse_index[search_name]):
                #  if successfully downloaded or already exists search for variant
                search_flag = self._unip.download_eve_data(name) if download_records else \
                    (search_name in self.eve_records)
                if search_flag:
                    res = self._eve_interperter(name, mutation, optimized)
                    if res != (-1, -1):
                        return res
                else:
                    continue
        return -1, -1

    def _eve_interperter(self, prot_name, mutation, optimized):
        """
        helper function for score_mutation_evemodel searches an evemodel csv for the mutation
        :param uniport_id: str Uid of the evemodel data
        :param mutation: Mutation obj to find
        :return: tuple (eve_score, -1 if not found, Eve_class 75% ASM prediction)
        """
        path = pjoin(EVE_VARIANTS_PATH, f"{prot_name}_HUMAN.csv")
        if not os.path.exists(path) or mutation.changeAA not in VALID_AA:
            return -1, -1
        if optimized == 1:
            data = load_parquet_cached(path, 'EVE')
        else:
            data = pd.read_csv(path, low_memory=False).fillna(-1)
        #  direct search first
        idx_predictions = data[data[EVE_POSITION_COLUMN] == mutation.loc]
        if len(idx_predictions) > 0 and idx_predictions[EVE_WT_COLUMN].unique()[0] == mutation.wt_aa:
            eve_variant = idx_predictions[idx_predictions[EVE_MUTATION_COLUMN] == mutation.changeAA]
            if pd.notna(eve_variant[EVE_SCORE_COLUMN].iloc[0]):
                return eve_variant.iloc[0][EVE_SCORE_COLUMN], eve_variant.iloc[0][EVE_PREDICTION_COLUMN]
        #  reference search
        references = mutation.ref_seqs.values()  # usually will be of size 1
        for reference in references:
            if len(reference) < REF_SEQ_PADDING * 2:  # skip if red seq is too short
                continue
            data_seq = ''.join(data[EVE_WT_COLUMN][0::20])
            seq_start = seacrh_by_ref_seq(data_seq, reference)
            if seq_start == -1:
                continue
            #  ref seq found in data - search through correct substitution
            all_substitution = data[seq_start:seq_start + N_AA]
            eve_variant = all_substitution[all_substitution[EVE_MUTATION_COLUMN] == mutation.changeAA]
            if len(eve_variant) != 1:
                warn_if(2, VERBOSE['thread_warnings'],
                        f"error with EVE search by reference, in {mutation.long_name} skipping...")
                continue
            else:
                if pd.notna(eve_variant[EVE_SCORE_COLUMN].iloc[0]):
                    return eve_variant.iloc[0][EVE_SCORE_COLUMN], eve_variant.iloc[0][EVE_PREDICTION_COLUMN]
        return -1, -1

    def score_mutation_eve_impute(self, mut, offset=0, gz=True, optimized=0):
        """
        :param mut: Mutation object
        :return: float CPT imputed score -1 if not found
        """
        entry_name = name_for_cpt(mut.protein.name)
        if entry_name == '': return -1, 'not_found'
        #  case 1 reference name found in CPT records
        if entry_name in self._cpt_ingene:
            score = self._eve_cpt_interperter(mut, entry_name, offset=offset, ingene=True, gz=gz)
            if score != -1: return score, 'cpt_ingene'
        if entry_name in self._cpt_exgene:
            score = self._eve_cpt_interperter(mut, entry_name, offset=offset, ingene=False, gz=gz)
            if score != -1: return score, 'cpt_exgene'
        #  case 2 expend search to all known protein references
        for name in mut.protein.aliases:
            entry_name = name_for_cpt(name)
            if entry_name in self._cpt_ingene:
                score = self._eve_cpt_interperter(mut, entry_name, offset=offset, ingene=True, gz=gz)
                if score != -1: return score, 'cpt_ingene'
            if entry_name in self._cpt_exgene:
                score = self._eve_cpt_interperter(mut, entry_name, offset=offset, ingene=False, gz=gz)
                if score != -1: return score, 'cpt_exgene'
        return -1, 'not_found'

    def _eve_cpt_interperter(self, mut, entry_name, ingene, offset, gz=True, optimized=0):
        """
        :param mut: Mutation object
        :param entry_name: Uniprot entery name
        :param ingene: bool search ingene or exgene records
        :param offset: int offset for mutation index
        :return: float score -1 if not found
        """
        desc = f"{mut.origAA}{mut.loc - offset}{mut.changeAA}"
        file = f"{entry_name}.csv"
        if ingene:
            path = pjoin(CPT_INGENE_PATH, file)
        else:
            path = pjoin(CPT_EXGENE_PATH, file)
        if gz:
            path += ".gz"
            df = pd.read_csv(path, compression='gzip', header=0, sep=',', quotechar='"', on_bad_lines='skip')
        else:
            df = pd.read_csv(path)

        direct_search = df[df[CPT_MUTATION_COLUMN] == desc][CPT_SCORE_COLUMN]
        if len(direct_search) > 0 and pd.notna(direct_search.iloc[0]):
            return float(direct_search.iloc[0])
        # mutation not found use reference sequence
        else:
            cpt_seq = sequence_from_cpt_df(df)
            references = mut.ref_seqs.values()  # usually will be of size 1
            for ref_seq in references:
                if len(ref_seq) < REF_SEQ_PADDING * 2:  # skip if ref seq is too short
                    continue
                row_idx = seacrh_by_ref_seq(main_seq=cpt_seq, ref_seq=ref_seq)
                if row_idx == -1:
                    continue  # ref not found
                else:
                    #  ref seq found in data - search through correct substitution
                    all_substitution = df[row_idx:row_idx + N_AA]
                    cpt_variant = all_substitution[all_substitution[CPT_MUTATION_COLUMN].str.endswith(mut.changeAA)]
                    if len(cpt_variant) != 1:
                        warn_if(2, VERBOSE['thread_warnings'],
                                f"error with CPT search by reference, in {mut.long_name} skipping...")
                        continue
                    elif pd.notna(cpt_variant.iloc[0][CPT_SCORE_COLUMN]):
                        return float(cpt_variant.iloc[0][CPT_SCORE_COLUMN])
            return -1

    def score_mutation_esm1b_precomputed(self, mut, offset=0, optimized=0):
        """
        :param mut: Mutation object
        :param offset: int offset to mutation index
        :return: float ESM-1b score, None if not found
        """
        entry_name = mut.protein.name
        if not entry_name:
            return None, 'failed'
        entry_name = name_for_esm(entry_name)
        #  case 1 reference name found in ESM records
        if entry_name in self.esm_index:
            search_name = self.esm_index[entry_name]
            score, score_type = self._esm_interperter(mut, search_name, offset, optimized=optimized)
            if score is not None:
                return score, score_type
        #  case 2 expend search to all known protein references
        #  case 3 unable to find - use reference sequence in isoforms
        temp_score = None
        for entry_name in mut.protein.aliases:
            entry_name = name_for_esm(entry_name)
            if entry_name in self.esm_index:
                search_name = self.esm_index[entry_name]
                score, score_type = self._esm_interperter(mut, search_name, offset, use_ref_seq=True,
                                                          optimized=optimized)
                if score is not None:
                    if score_type == 'direct':
                        return score, score_type
                    else:
                        #  may override previous case of multiple temp_scores is rare and not trivial to handle
                        # always prefers direct scores
                        temp_score = score
        if temp_score:
            return temp_score, score_type
        return None, 'failed'

    @staticmethod
    def score_mutation_esm_inference(model, alphabet, mut, method='masked_marginals', offset=1, log='inference'):
        """
        :param model: esm initiated model
        :param alphabet: alphabet api
        :param mut: Mutation object
        :param offset: int offset to mutation index
        :param log: str
        :param method: scoring method str: wt_marginals | mutant_marginals | masked_marginals
        :return: tuple (float | None ESM3 masked marginal, str log)
        """
        if not mut.isoforms:
            return None, log + '\tno_iso'
        if mut.changeAA == 'X':  # invalid change
            return None, None
        mut_idx = mut.loc - offset
        log = log + '\toriginal_sequence'
        sequence = mut.sequence(how='min')
        if len(sequence) > ESM_MAX_LENGTH:
            new_offset, sequence = esm_process_long_sequences(sequence, mut_idx)
            mut_idx = mut_idx - new_offset
            log = log + '\ttrimmed_sequence'
        assert sequence[mut_idx] == mut.origAA, 'sequence does not match wild-type AA - check offset'
        if method == 'mutant_marginals':
            input_seq = sequence[:mut_idx] + mut.changeAA + sequence[mut_idx + 1:]
        elif method == 'masked_marginals':
            input_seq = sequence[:mut_idx] + MASK_TOKEN + sequence[mut_idx + 1:]
        else:
            input_seq = sequence
        tokenizer = alphabet.get_batch_converter()
        # tokens = torch.tensor(tokenized, dtype=torch.int64).unsqueeze(0).to(DEVICE)
        _, _, batch_tokens = tokenizer([(mut.long_name, input_seq)])
        batch_tokens = batch_tokens.to(DEVICE)
        logits = esm_seq_logits(model=model, tokens=batch_tokens, log=True, softmax=True, return_device='cpu',
                                esm_version=1)
        # use masked marginals score
        return float(logits[mut_idx][AA_ESM_LOC[mut.changeAA]] - logits[mut_idx][AA_ESM_LOC[mut.origAA]]), log

    @staticmethod
    def score_mutation_esm3(model, tokenizer, mut, method='masked_marginals', offset=1, log=''):
        """
        :param model: ESM3 initiated model
        :param tokenizer: ESM3 sequence tokenizer
        :param mut: Mutation object
        :param offset: int offset to mutation index
        :param log: str
        :param method: scoring method str: wt_marginals | mutant_marginals | masked_marginals
        :return: tuple (float | None ESM3 masked marginal, str log)
        """
        if not mut.isoforms:
            return None, log + '\tno_iso'
        mut_idx = mut.loc - offset
        log = log + '\toriginal_sequence'
        #  choose sequence of minimal length
        sequence = mut.sequence(how='min')
        if len(sequence) > ESM_MAX_LENGTH:
            new_offset, sequence = esm_process_long_sequences(sequence, mut_idx)
            mut_idx = mut_idx - new_offset
            log = log + '\ttrimmed_sequence'
        #  mask sequence
        assert sequence[mut_idx] == mut.origAA, 'sequence does not match wild-type AA - check offset'
        # wt_marginals is default
        input_seq = sequence
        if method == 'mutant_marginals':
            input_seq = sequence[:mut_idx] + mut.changeAA + sequence[mut_idx + 1:]
        elif method == 'masked_marginals':
            input_seq = sequence[:mut_idx] + '_' + sequence[mut_idx + 1:]
        tokenized = tokenizer(input_seq)['input_ids']
        tokens = torch.tensor(tokenized, dtype=torch.int64).unsqueeze(0).to(DEVICE)
        print(f"calculating score {mut.long_name}")
        logits = esm_seq_logits(model=model, tokens=tokens, log=True, softmax=True, return_device='cpu')
        # use masked marginals score
        return logits[mut_idx][AA_ESM_LOC[mut.changeAA]] - logits[mut_idx][AA_ESM_LOC[mut.origAA]], log

    @staticmethod
    def _esm_interperter(mut, search_name, offset, use_ref_seq=False, optimized=0):
        file_path = pjoin(ESM_VARIANTS_PATH, search_name + ESM_FILE_SUFFIX)
        column = f"{mut.origAA} {mut.loc - offset}"
        if optimized == 1:
            data = load_parquet_cached(file_path, 'ESM')
        else:
            data = pd.read_csv(file_path)
        assert mut.changeAA in A_A, f'invalid Amino Acid symbol {mut.long_name}'
        if mut.changeAA == 'X':  # invalid change
            return None, None
        row = data[data.iloc[:, 0] == mut.changeAA].index.tolist()[0]
        if column in data.columns:
            return float(data[column][row]), 'direct'
        if use_ref_seq:
            esm_seq = sequence_from_esm_df(data)
            references = mut.ref_seqs.values()
            for ref in references:
                seq_start = esm_seq.find(ref)
                if seq_start == -1:
                    continue
                # since seq ignores first column add 1
                aa_location = esm_seq.find(mut.origAA, seq_start, seq_start + 10) + 1
                column = data.iloc[:, aa_location]
                return float(column[row]), 'indirect'
        return None, None

    def score_mutation_afm(self, mutation, chunk=None, uid_index=None, explicitly_access=None, offset=0,
                           use_alias=False, use_unreviewed_uids=False):
        """
        Scores mutation using AlphaMissense
        Searches main Uniprot entry and if not found all reviewed entries

        :param mutation: Mutation object
        :param offset: int offset for mutation index
        :param chunk: DataFrame - optional if given will search for the mutation in batch
                    instead of importing data from disk. Recommended when scoring multiple mutations.
        :param uid_index: optional set of uids in chunk, can reduce running time
        :param explicitly_access: optional set of uids that should explicitly imported from the DB
        :param use_alias: optional bool if true will use alias uids in search
        :param use_unreviewed_uids: optional bool if true will use unreviewed alias uids
        :return: float AlphaMissense score if found else -1
        """
        main_uid = mutation.protein.Uid
        alias_uids = mutation.protein.all_uids()['reviewed'] if not use_unreviewed_uids else \
            mutation.protein.all_uids()['non_reviewed']
        references = mutation.ref_seqs.values()
        if isinstance(alias_uids, str):
            alias_uids = [alias_uids]
        variant_location = mutation.loc + offset
        variant = f"{mutation.origAA}{variant_location}{mutation.changeAA}"
        index = uid_index if uid_index else self.af_index
        explicitly_access = explicitly_access if explicitly_access else set()
        uids = [main_uid] + alias_uids if use_alias else [main_uid]
        main_uid_flag = True
        for uid in uids:
            if uid in index:
                #  search by chunks will be removed in the future
                if chunk is not None:
                    score, desc = self._chunk_score_afm(chunk, uid, variant, explicitly_access, references)
                else:  # default search by AFM protein files
                    score, desc = self._afm_score_from_uid(uid, variant, mut_references=references)

                if score is not None and main_uid_flag:
                    return score, f'main_uid\t{desc}'
                elif score is not None and not main_uid_flag:
                    return score, f'alias\t{desc}'
                main_uid_flag = False
        # score not found
        return None, NO_TYPE

    #  Will be removed inb the future
    def _chunk_score_afm(self, chunk, uid, variant, explicitly_access, references):
        idx_from = self.af_index[uid]
        idx_to = self.af_ranges[str(idx_from)]
        if uid in explicitly_access:
            # used for sequences that spread over multiple chunks --> load ranges explicitly
            score, desc = self._afm_score_from_uid(uid, variant, idx_from, idx_to,
                                                   mut_references=references, chunk=None)
        else:
            score, desc = self._afm_score_from_uid(uid, variant, idx_from, idx_to,
                                                   mut_references=references, chunk=chunk)
        return score, desc

    def _afm_score_from_uid(self, uid, variant, idx_from=None, idx_to=None, mut_references=None, chunk=None):
        """
        extract alpha Missense score from raw data if not found returns -1
        :param variant: str protein variant in formal [AA][location][AA]
        :param idx_from: int
        :param idx_to: int
        :param v: DataFrame - optional if given will search for the mutation in batch
                instead of importing data from disk. Recommended when scoring multiple mutations.
        :return:
        """
        if not (variant[-1] in VALID_AA) and (variant[0] in VALID_AA):
            return None, ''
        if chunk is not None:
            data = chunk
        elif idx_to is not None and idx_from is not None:
            data = afm_range_read(idx_from, idx_to)
        else:
            path = pjoin(AFM_PROTEINS_PATH, f'{uid}.csv.gz')
            if not os.path.exists(path):
                return None, ''
            data = pd.read_csv(path, compression='gzip', header=0, sep=',', quotechar='"', on_bad_lines='skip')
        score = data[(data['protein_variant'] == variant) & (data['uniprot_id'] == uid)]['am_pathogenicity']
        if not score.empty:
            return float(score.iloc[0]), 'direct'
        # direct search failed --> try by reference
        if mut_references:
            record_data = data[data['uniprot_id'] == uid]
            afm_seq = sequence_from_afm_df(record_data)
            if not afm_seq:
                warn_if(2, VERBOSE['thread_warnings'], f'afm seq truncated {uid}')
                record_data.to_csv(f'log_{variant}_{uid}.csv')
                return None, ''
            for ref_seq in mut_references:
                if len(ref_seq) < REF_SEQ_PADDING * 2:  # skip if ref seq is too short
                    continue
                row_idx = seacrh_by_ref_seq(main_seq=afm_seq, ref_seq=ref_seq, n_rows=19)
                if row_idx == -1:
                    continue  # ref not found
                else:
                    #  ref seq found in data - search through correct substitution
                    all_substitution = record_data[row_idx:row_idx + 19]
                    afm_variant = all_substitution[all_substitution['protein_variant'].str.endswith(variant[-1])]
                    if len(afm_variant) != 1:
                        warn_if(2, VERBOSE['thread_warnings'],
                                f"error with AFM search by reference, in {variant} skipping...")
                        record_data.to_csv(f'log_{variant}_{uid}.csv')
                        continue
                    elif pd.notna(afm_variant.iloc[0]['am_pathogenicity']):
                        return float(afm_variant.iloc[0]['am_pathogenicity']), 'ref'
        return None, ''

    def score_interface(self, mutation, use_alias=False, use_unreviewed_uids=False):
        """
        interface score using pioneer
        :param mutation: Mutation object
        :return: tuple (int score, str score_type)
        """
        if self.pioneer_predictions is None or self.pioneer_sequences is None:
            self._load_pioneer_data()
        main_uid = mutation.protein.Uid
        alias_uids = mutation.protein.all_uids()['reviewed'] if not use_unreviewed_uids else \
            mutation.protein.all_uids()['non_reviewed']
        references = mutation.ref_seqs.values()
        if isinstance(alias_uids, str):
            alias_uids = [alias_uids]
        uids = [main_uid] + alias_uids if use_alias else [main_uid]
        for uid in uids:
            if uid not in self.pioneer_sequences:
                continue
            pioneer_seqs = self.pioneer_sequences[uid]
            #  TODO search only by ref does not cover edges mutation.loc < 10 or mutation.loc > len(seq)-10
            for ref_seq in references:
                if len(ref_seq) < REF_SEQ_PADDING * 2:  # skip if ref seq is too short
                    continue
                for pioneer_seq in pioneer_seqs:  # usually of size 1
                    pred_idx = seacrh_by_ref_seq(pioneer_seq, ref_seq, n_rows=1)
                    if pred_idx == -1:  # ref not found in pioneer seq
                        continue
                    pred = self.pioneer_predictions[pioneer_seq][pred_idx]
                    return pred
        return None

