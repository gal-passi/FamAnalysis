import os
from os.path import basename
import subprocess
from subprocess import check_output
import numpy as np
from Bio import SeqIO
import json
import warnings
import glob
import pandas as pd
import Mutation
import pickle
from utils import adaptive_chunksize, make_fasta, afm_range_read
import Connections
from pathlib import Path
import Patient
import Family
from copy import deepcopy
from definitions import *
import deprecation

ROOT = os.path.join(os.path.dirname(__file__))


class ProteinAnalyzer:
    INTERFACE = "/cs/staff/dina/utils/srcs/interface/interface"
    STRINGS_REVERSE_INDEX_ROOT = f"{ROOT}/DB/strings/strings_rindex.txt"
    PDBS_DB_ROOT = f"{ROOT}/DB/pdbs"

    def __init__(self, *proteins, verbose_level=1):
        self.proteins = {protein.name: protein for protein in proteins}
        self._cpt_ingene = {basename(p)[:-CPT_TRUNCATE_TAIL] for p in glob.glob(pjoin(CPT_INGENE_PATH, '*.csv.gz'))}
        self._cpt_exgene = {basename(p)[:-CPT_TRUNCATE_TAIL] for p in glob.glob(pjoin(CPT_EXGENE_PATH, '*.csv.gz'))}
        self.eve_records = {basename(p)[:-EVE_TRUNCATE_TAIL] for p in glob.glob(pjoin(EVE_VARIANTS_PATH, '*.csv'))}
        self._unip = Connections.Uniport(verbose_level=verbose_level)
        with open(self.STRINGS_REVERSE_INDEX_ROOT, "rb") as fp:
            self.strings_rindex = pickle.load(fp)
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

    def analyze_all_proteins(self):
        """
        analyzes all protein's mutations
        :return: {name_mutation: Scores obj}
        """
        ret = {}
        for protein in self.proteins:
            ret = ret | self.analyze_single_protein(self.proteins[protein])

        return ret

    def analyze_single_protein(self, protein):
        """
        analyzes a single Protein obj
        :param protein: Prot obj
        :return: {name_mutation: Scores obj}
        """
        ret = {}
        for mut in protein.muts:
            ret = ret | self.analyze_single_mutation(mut, protein)

        return ret

    def analyze_single_mutation(self, mutation, protein):
        """
        analyzes a single mutation of a Protein obj
        :param mutation: p.{AA}{location}{AA} must be in prteins mutation dictionary
        :param protein: Protein obj.
        :return: {name_mutation_isoform: Scores obj}
        """
        paths = protein.pdb_paths(mutation)  # {isoform: [pdb_pathes]}
        extended_name = False
        if len(paths) > 1:
            warnings.warn("Notice: mutation found in more then two isoforms, will be saved in format name_mut_iso")
            extended_name = True
        for iso, pdbs in paths:
            for pdb in pdbs:
                chains = self.pdb_chains(pdb)  # {chain id: Bio.seq}
                isoform = iso if extended_name else ""
                ref = protein.get_ref_seq(mutation, isoform)
                relevant_ids = self._find_chains(chains, ref)  # set
                self.score_distance(pdb, chains, relevant_ids, mutation)

    def score_distances(self, mutation):
        """
        :param mutation: Mutation object
        :return: {pdb: interface_distance}
        """
        scores = {}
        for pdbid in mutation.raw_pdbs:
            path = os.path.join(self.PDBS_DB_ROOT, f"pdb{pdbid}.ent")
            chains = self.pdb_chains(path)
            mutation_chains = self._find_chains(chains, mutation)
            scores[pdbid] = self._score_distance(path, set(chains.keys()), mutation_chains)
        return scores

    #  TODO finish - we should preforme the check to all pairs of relevent_id - other chains in diffiernt thr
    #  TODO start form the lowest thr and move up - if lowest was match no reason to continue
    #  TODO no reason the check pairs from two different directions AC == CA ids are a set use it!
    def _score_distance(self, pdb_path, all_chains, mutation_chains):
        """
        gives a distance score to the given mutation in the pdb_path.
        :param pdb_path: path to pdb file in question
        :param chains: {chain_id, seq}
        :param ids: relevent ids for mutation
        :param ids: threshold int threshold parameter for interface script
        :param mutation: mutation object
        :return: int 0-3 where 0 is very far and 3 is very close
        """
        interface_score = np.inf
        for d in ['2', '3', '4', '6', '8', '12', '16']:
            for chain_id, loc in mutation_chains:
                search_term = b'%d  %b' % (loc, chain_id.encode('utf8'))
                other_chains = all_chains.difference(chain_id)
                for c in other_chains:
                    res = check_output([self.INTERFACE, '-c', pdb_path, chain_id, c, d]).find(search_term)
                    if res != -1:
                        return d
        return interface_score

        # start from location 2

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
    def score_mutation_evemodel(self, protein, mutation, prot_name=None, download_records=False):
        """
        looks for an evemodel score for the mutation given
        :param download_records: bool whether to download new eve records. This should be True is user choose setup
               with partial eve data.
        :param protein: Protein object
        :param mutation: Mutation object
        :return: (score, prediction) if not found (-1,-1)
        """
        search_name = protein.name if not prot_name else prot_name
        #  Case 1: Protein reference in is directly found in eve records
        if search_name in self.eve_records:
            res = self._eve_interperter(search_name, mutation)
            if res != (-1, -1):
                return res
        # Case 2: Protein has an alias directly found in eve records
        for alias in protein.aliases:
            if alias in self.eve_records:
                res = self._eve_interperter(alias, mutation)
                if res != (-1, -1):
                    return res
        #  Case 3: Protein main reference is an alias of a record directly covered by eve
        if search_name in self.eve_reverse_index:
            for name in set(self.eve_reverse_index[search_name]):
                #  if successfully downloaded or already exists search for variant
                search_flag = self._unip.download_eve_data(name) if download_records else \
                    (search_name in self.eve_records)
                if search_flag:
                    res = self._eve_interperter(name, mutation)
                    if res != (-1, -1):
                        return res
                else:
                    continue
        return -1, -1

    def _eve_interperter(self, prot_name, mutation):
        """
        helper function for score_mutation_evemodel searches an evemodel csv for the mutation
        :param uniport_id: str Uid of the evemodel data
        :param mutation: Mutation obj to find
        :return: tuple (eve_score, -1 if not found, Eve_class 75% ASM prediction)
        """
        path = pjoin(EVE_VARIANTS_PATH, f"{prot_name}_HUMAN.csv")
        if not os.path.exists(path):
            return -1, -1
        data = pd.read_csv(path, low_memory=False).fillna(-1)
        references = mutation.ref_seqs.values()  # usually will be of size 1
        for reference in references:
            data_seq = ''.join(data["wt_aa"][0::20])
            seq_start = data_seq.find(reference)  # TODO handle edges - will be implemented in Mutation
            if seq_start == -1:
                continue
            aa_location = data_seq.find(mutation.origAA, seq_start, seq_start + 10)
            row_index = (aa_location * 20) + AA_TO_INDEX_EVE[mutation.changeAA]
            return data.iloc[row_index]["EVE_scores_ASM"], data.iloc[row_index]["EVE_classes_75_pct_retained_ASM"]
        return -1, -1

    def score_mutation_eve_impute(self, mut, offset=0, gz=True):
        """
        :param mut: Mutation object
        :return: float CPT imputed score -1 if not found
        """
        name = f"{mut.origAA}{mut.loc - offset}{mut.changeAA}"
        entery_name = mut.protein.entery_name()
        if entery_name == '': return -1, 'not_found'
        #  case 1 reference name found in CPT records
        if entery_name in self._cpt_ingene:
            score = self._eve_cpt_interperter(name, entery_name, ingene=True, gz=gz)
            if score != -1: return score, 'cpt_ingene'
        if entery_name in self._cpt_exgene:
            score = self._eve_cpt_interperter(name, entery_name, ingene=False, gz=gz)
            if score != -1: return score, 'cpt_exgene'
        #  case 2 expend search to all known protein references
        for entery_name in mut.protein.entery_name(all=True):
            if entery_name in self._cpt_ingene:
                score = self._eve_cpt_interperter(name, entery_name, ingene=True, gz=gz)
                if score != -1: return score, 'cpt_ingene'
            if entery_name in self._cpt_exgene:
                score = self._eve_cpt_interperter(name, entery_name, ingene=False, gz=gz)
                if score != -1: return score, 'cpt_exgene'
        return -1, 'not_found'

    def _eve_cpt_interperter(self, desc, entery_name, ingene, gz=True):
        """
        :param desc: missmatch description in form {AA}{index}{AA}
        :param entery_name: Uniprot entery name
        :param ingene: bool search ingene or exgene records
        :return: float score -1 if not found
        """
        # print(f"called with {entery_name} and ingene={ingene}")
        file = f"{entery_name}.csv"
        if ingene:
            path = pjoin(CPT_INGENE_PATH, file)
        else:
            path = pjoin(CPT_EXGENE_PATH, file)
        if gz:
            path += ".gz"
            df = pd.read_csv(path, compression='gzip', header=0, sep=',', quotechar='"', on_bad_lines='skip')
        else:
            df = pd.read_csv(path)

        df = df[df[CPT_MUTATION_COLUMN] == desc][CPT_SCORE_COLUMN]
        return -1 if len(df) == 0 else float(df.iloc[0])

    def score_mutation_afm(self, mutation, chunk=None, uid_index=None, offset=0, use_alias=False):
        """
        tries to give scores for mutation using AlphaMissense
        will search main Uniprot entry and if not found all reviewed entries

        :param mutation: Mutation object
        :param offset: int offset for mutation index
        :param chunk: DataFrame - optional if given will search for the mutation in batch
                    instead of importing data from disk. Recommended when scoring multiple mutations.
        :param uid_index: optional set of uids in chunk, can reduce running time
        :return: float AlphaMissense score if found else -1
        """
        main_uid = mutation.protein.Uid
        reviewed_uids = mutation.protein.all_uids()['reviewed']
        if isinstance(reviewed_uids, str):
            reviewed_uids = [reviewed_uids]
        variant_location = mutation.loc + offset
        variant = f"{mutation.origAA}{variant_location}{mutation.changeAA}"
        index = uid_index if uid_index is not None else self.af_index
        uids = [main_uid] + reviewed_uids if use_alias else [main_uid]
        for uid in uids:
            if uid in index:
                idx_from = self.af_index[uid]
                idx_to = self.af_ranges[str(idx_from)]
                score = self._afm_score_from_uid(variant, idx_from, idx_to, chunk)
                if score is not None:
                    return score
        # score not found
        return None

    def score_mutation_esm(self, mut, offset=0):
        """
        :param mut: Mutation object
        :return: float ESM-1b score, None if not found
        """
        entery_name = mut.protein.entery_name()
        if entery_name == '':
            return None
        #  case 1 reference name found in ESM records
        if entery_name in self.esm_index:
            search_name = self.esm_index[entery_name]
            score = self._esm_interperter(mut, search_name, offset)
            if score is not None:
                return score
        #  case 2 expend search to all known protein references
        for entery_name in mut.protein.entery_name(all=True):
            if entery_name in self.esm_index:
                search_name = self.esm_index[entery_name]
                score = self._esm_interperter(mut, search_name, offset)
                if score is not None:
                    return score
        return None

    @staticmethod
    def _esm_interperter(mut, search_name, offset):
        file_path = pjoin(ESM_VARIANTS_PATH, search_name + ESM_FILE_SUFFIX)
        row = AA_TO_INDEX_ESM[mut.changeAA]
        column = f"{mut.origAA} {mut.loc - offset}"
        data = pd.read_csv(file_path)
        if column in data.columns:
            return float(data[column][row])
        else:
            return None

    @staticmethod
    def _afm_score_from_uid(variant, idx_from, idx_to, chunk=None):
        """
        extract alpha Missense score from raw data if not found returns -1
        :param variant: str protein variant in formal [AA][location][AA]
        :param idx_from: int
        :param idx_to: int
        :param v: DataFrame - optional if given will search for the mutation in batch
                    instead of importing data from disk. Recommended when scoring multiple mutations.
        :return:
        """
        data = afm_range_read(idx_from, idx_to) if chunk is None else chunk
        score = data[data['protein_variant'] == variant]['am_pathogenicity']
        return float(score.iloc[0]) if not score.empty else None

    @staticmethod
    def mutations_by_rank():
        """
        :return: 4 sets of mutations (r1, r2, r3,r4)
        """
        n, c = len(glob.glob('DB/Patients/HR*')), 0
        r1, r2, r3, r4 = set(), set(), set(), set()
        for p in [os.path.basename(i)[:-4] for i in glob.glob('DB/Patients/HR*')]:
            patient = Patient.Patient(p)
            f_id = p.split('_')[0]
            family = Family.Family(family_id=f_id)
            rank_1 = patient.top_all() | set(patient.in_literature())
            recurring = family.intersect_all()
            rank_2 = recurring
            rank_3 = patient.in_n_params(2) | set(patient.top_interface(in_literature=True))
            rank_4 = (set(patient.top_interface(in_literature=True)) & patient.in_n_params(1, interface=0)) | \
                     patient.in_n_params(3) | (set(patient.in_literature()) & patient.in_n_params(2))
            if recurring:
                m_counter = {m: [] for m in recurring}
                for f in [os.path.basename(i)[:-4] for i in glob.glob('DB/Family/*')]:
                    temp_family = Family.Family(family_id=f)
                    common_muts = recurring & temp_family.intersect_all()
                    for m in common_muts:
                        m_counter[m].append(f)
                with open(f'family_{f_id}_recurring.txt', 'w+') as convert_file:
                    write_form = {k.long_name: v for k, v in m_counter.items()}
                    convert_file.write(json.dumps(write_form))
                m_counter = [k for k, v in m_counter.items() if len(v) > 1]
                rank_4 = rank_4 | set(m_counter)

            all_muts = list(deepcopy(rank_1))
            rank_1 = rank_1 - (rank_2 | rank_3 | rank_4)
            rank_2 = rank_2 - (rank_3 | rank_4)
            rank_3 = rank_3 - rank_4
            r1 = r1 | rank_1
            r2 = r2 | rank_2
            r3 = r3 | rank_3
            r4 = r4 | rank_4
            c += 1
            print(f"finished {round(c / n, 2) * 100}%")
        print("done")
        r1 = r1 - (r2 | r3 | r4)
        r2 = r2 - (r3 | r4)
        r3 = r3 - r4
        return r1, r2, r3, r4
