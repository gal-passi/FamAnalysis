import os
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
    BERT_DB_ROOT = f"{ROOT}/DB/bert"
    STRINGS_DB_ROOT = f"{ROOT}/DB/strings"
    STRINGS_INDEX_ROOT = f"{ROOT}/DB/strings/strings_index.txt"
    STRINGS_REVERSE_INDEX_ROOT = f"{ROOT}/DB/strings/strings_rindex.txt"
    STRING_NEBS = f"{ROOT}/DB/strings/nebs.txt"
    ESM_MODEL = "/cs/labs/dina/gal_passi/breast_cancer/pdbs-venv/esm"
    ALPHAFOLD_PREDICTIONS = f"{ROOT}/DB/alphfold_out"
    ALPHAFOLD_JOBS = f"{ROOT}/temp_alphfold_jobs"
    MUTATION_PAIRS = f"{ROOT}/"
    PDBS_DB_ROOT = f"{ROOT}/DB/pdbs"
    AA_TO_INDEX = {"A": 0, "C": 1, "D": 2, "E": 3, "F": 4, "G": 5, "H": 6, "I": 7, "K": 8, "L": 9, "M": 10,
                   "N": 11, "P": 12, "Q": 13, "R": 14, "S": 15, "T": 16, "V": 17, "W": 18, "Y": 19}

    def __init__(self, *proteins, verbose_level=1):
        self.proteins = {protein.name: protein for protein in proteins}
        self._cpt_ingene = set([Path(i).with_suffix('').stem for i in glob.glob(f"{pjoin(CPT_INGENE_PATH, '*.csv*')}")])
        self._cpt_exgene = set([Path(i).with_suffix('').stem for i in glob.glob(f"{pjoin(CPT_EXGENE_PATH, '*.csv*')}")])
        self.eve_records = list(map(lambda x: os.path.basename(x)[:-10], glob.glob(EVE_PATH + "*.csv")))
        #with open(EVE_INDEX_PATH, "rb") as fp:   # Unpickling
        #    self.eve_index_unip = pickle.load(fp)
        #with open(EVE_EXTENDED_INDEX_PATH, "rb") as fp:   # Unpickling
        #    self.eve_index_unip_full = pickle.load(fp)
        with open(self.STRINGS_INDEX_ROOT, "rb") as fp:
            self.strings_index = pickle.load(fp)
        with open(self.STRINGS_REVERSE_INDEX_ROOT, "rb") as fp:
            self.strings_rindex = pickle.load(fp)
        with open(self.STRING_NEBS, "rb") as fp:
            self.protein_nebs = pickle.load(fp)
        with open(EVE_INDEX_PATH_2, "rb") as fp:
            self.eve_index = pickle.load(fp)
        with open(EVE_INVERSE_INDEX, "rb") as fp:
            self.eve_reverse_index = pickle.load(fp)
        self._unip = Connections.Uniport(verbose_level= verbose_level)
        with open(AFM_DIRECTORY_PATH, 'r') as file:
            self.af_index = json.load(file)
        with open(AFM_RANGES_PATH, 'r') as file:
            self.af_ranges = json.load(file)


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
        for d in ['2','3','4','6','8','12','16']:
            for chain_id, loc in mutation_chains:
                search_term = b'%d  %b'%(loc, chain_id.encode('utf8'))
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
            res = {(c_id, seq.find(ref)+offset) for c_id, seq in chains_dict.items() if ref in seq}
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
        return {id.split(':')[1] if ':' in id else id : seq for id, seq in res.items()}

    @staticmethod
    def search_eve_record(prot_name):
        """
        searches for an existing evemodel data
        :param prot_name as in Protein.name:
        :return: bool
        """
        return prot_name in list(map(lambda x: os.path.basename(x)[:-10], glob.glob("DB/EVE/variant_files/")))

    # TODO change eve db to work with unip - search by name might miss some proteins
    def score_mutation_evemodel(self, protein, mutation, prot_name=None):
        """
        looks for an evemodel score for the mutation given
        :param protein: Protein object
        :param mutation: Mutation object
        :return: (score, prediction) if not found (-1,-1)
        """
        search_name = protein.name if not prot_name else prot_name
        #print(f"searching using {search_name}")
        if search_name in self.eve_records:
            res = self._eve_interperter(search_name, mutation)
            if res != (-1, -1):
                return res
        for alias in self._unip.synonms(protein):
            if alias in self.eve_records:
                res = self._eve_interperter(alias, mutation)
                if res != (-1, -1):
                    return res
        if search_name in self.eve_reverse_index:
            for eve_data in set(self.eve_reverse_index[search_name]):
                #print(f"trying with new searching using {eve_data}")
                res = self._eve_interperter(eve_data, mutation)
                if res != (-1, -1):
                    return res
        return -1, -1

    def _eve_interperter(self, prot_name, mutation):
        """
        helper function for score_mutation_evemodel searches an evemodel csv for the mutation
        :param uniport_id: str Uid of the evemodel data
        :param mutation: Mutation obj to find
        :return: tuple (eve_score, -1 if not found, Eve_class 75% ASM prediction)
        """
        path = EVE_PATH + f"{prot_name}_HUMAN.csv"
        if not os.path.exists(path):
            return -1, -1
        data = pd.read_csv(path).fillna(-1)
        references = mutation.ref_seqs.values()  # usually will be of size 1
        for reference in references:
            data_seq = ''.join(data["wt_aa"][0::20])
            #print(data_seq)
            seq_start = data_seq.find(reference)  # TODO handle edges - will be implemented in Mutation
            if seq_start == -1:
                continue
            aa_location = data_seq.find(mutation.origAA, seq_start, seq_start + 10)
            row_index = (aa_location * 20) + self.AA_TO_INDEX[mutation.changeAA]
            # For debugging purposes
            # print(row_index)
            # print(data.iloc[row_index]["mt_aa"])
            return data.iloc[row_index]["EVE_scores_ASM"], data.iloc[row_index]["EVE_classes_75_pct_retained_ASM"]
        return -1, -1

    def Bert_score(self, prot, path=BERT_DB_ROOT, mut_name="", isoform="", seq="", log={"success": [], "fail_iso": [], "no_mut": []}, offset=0):
        """
        calculates Bert scors using models esm1v_t33_650M_UR90S  and saves to csv at path
        this might take up to 15 minutes per mutation
        :param prot: Protein obj
        :param mut_name: Optional str to run on one mutation at a time, should be same as in protein.mutations keys
        :param isoform: Optional specific isoform to run on
        """

        def _handle_bert(prot_obj, m_name, final_seq, bert_offset, csv_path):
            path = os.path.join(csv_path, prot_obj.name +'_'+ m_name + ".csv")
            df = pd.DataFrame({'mutation': [Mutation.Mutation(m_name, prot_obj, prot_obj.mutations[m_name]).for_bert_location]})
            df.to_csv(path)
            script_location = os.path.join(self.ESM_MODEL, "variant-prediction/predict.py")

            # NOTE CHANGE script_location IF RUN FROM DIFFERENT LOCATION
            print("Activating BERT")
            with subprocess.Popen(f"python3 {script_location} "
                                  f"--model-location esm1v_t33_650M_UR90S_1 esm1v_t33_650M_UR90S_2 esm1v_t33_650M_UR90S_3 esm1v_t33_650M_UR90S_4 esm1v_t33_650M_UR90S_5 "
                                  f"--sequence {final_seq} "
                                  f"--dms-input {path} "
                                  f"--mutation-col mutation "
                                  f"--dms-output {path} "
                                  f"--offset-idx {bert_offset} "
                                  f"--scoring-strategy wt-marginals", stdout=subprocess.PIPE, shell=True) as proc:
                print(proc.stdout.read())
                print(f"done {prot_obj.name}({m_name}) - len = {len(final_seq)}, off = {bert_offset}")
                log['success'].append(prot.name + "_" + m_name)
                return log

        if not prot.mutations:
            print(f"done {prot.name} no mutations found to process")
            log['no_mut'].append(prot.name)
            return log
        if not seq:
            for name, data in prot.mutations.items():
                mut = Mutation.Mutation(name, prot, data)
                if len(mut.isoforms) == 0:
                    warnings.warn(f"No valid isoforms for {prot.name}({mut.name})")
                    log['fail_iso'].append(prot.name)
                    continue
                iso = next(iter(mut.isoforms)) if not isoform else isoform
                if (len(mut.isoforms) > 1) and (not isoform):
                    warnings.warn(f"More then 1 isoform found using {iso} use isoform parameter to choose a different one from {mut.isoforms}")
                sequence = prot.isoforms[iso]

                if len(sequence) > 1020:
                    offset, sequence = self._bert_process_long_sequences(sequence, mut)

                log = _handle_bert(prot, name, sequence, offset, path)

        if seq and mut_name:
            sequence = seq
            return _handle_bert(prot, mut_name, sequence, offset, path)

    @staticmethod
    def _bert_process_long_sequences(seq, mut):
        """
        trims seq to len < 1024 using mut Object. to fit for bert model
        """
        loc = mut.loc
        left_bound = 0 if loc-510 < 0 else loc - 510
        right_bound = len(seq) if loc + 510 >len(seq) else loc+510
        left_excess = 0 if loc - 510 > 0 else abs(loc-510)
        right_excess = 0 if loc + 510 <= len(seq) else loc+510-len(seq)
        if (left_excess == 0) and (right_excess == 0):
            #print("option 1")
            return left_bound, seq[left_bound:right_bound]
        if left_excess > 0:
            #print("option 2")
            return 0, seq[left_bound:right_bound+left_excess]
        if right_excess > 0:
            #print("option 3")
            return left_bound-right_excess, seq[left_bound-right_excess:right_bound]

    def score_mutation_eve_impute(self, mut, offset = 0, gz=False):
        """
        :param mut: Mutation object
        :return: float CPT imputed score -1 if not found
        """
        name = f"{mut.origAA}{mut.loc - offset}{mut.changeAA}"
        entery_name = mut.protein.entery_name()
        if entery_name == '': return -1
        if entery_name in self._cpt_ingene:
            #print("found in ingene")
            score = self._eve_cpt_interperter(name, entery_name, ingene=True, gz=gz)
            if score != -1: return score
        if entery_name in self._cpt_exgene:
            #print("found in ingene")
            score = self._eve_cpt_interperter(name, entery_name, ingene=False, gz=gz)
            if score != -1: return score
        for entery_name in mut.protein.entery_name(all=True):
            if entery_name in self._cpt_ingene:
                #print(f"found {entery_name} in ingene (2)")
                score = self._eve_cpt_interperter(name, entery_name, ingene=True, gz=gz)
                if score != -1: return score
            if entery_name in self._cpt_exgene:
                #print(f"found {entery_name} in exgene (2)")
                score = self._eve_cpt_interperter(name, entery_name, ingene=False, gz=gz)
                if score != -1: return score
        return -1

    def _eve_cpt_interperter(self, desc, entery_name, ingene, gz=False):
        """
        :param desc: missmatch description in form {AA}{index}{AA}
        :param entery_name: Uniprot entery name
        :param ingene: bool search ingene or exgene records
        :return: float score -1 if not found
        """
        #print(f"called with {entery_name} and ingene={ingene}")
        file = f"{entery_name}.csv"
        if ingene:
            path = pjoin(CPT_INGENE_PATH, file)
        else:
            path = pjoin(CPT_INGENE_PATH, file)
        if gz:
            path += ".gz"
            df = pd.read_csv(path, compression='gzip', header=0, sep=',', quotechar='"', on_bad_lines='skip')
        else:
            df = pd.read_csv(path)

        df = df[df['mutant'] == desc][f'{"Transfer_EVE" if ingene else "Transfer_imputed_EVE"}']
        return -1 if len(df) == 0 else float(df)

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
            r2 = r2 |rank_2
            r3 = r3 | rank_3
            r4 = r4 | rank_4
            c+=1
            print(f"finished {round(c / n, 2) * 100}%")
        print("done")
        r1 = r1 - (r2 | r3 | r4)
        r2 = r2 - (r3 | r4)
        r3 = r3 - r4
        return r1, r2, r3, r4

    def find_mutation_pairs(self, patient):  #TODO not yet implemented
        """
        pairs of interacting proteins which are both damaged in the same patient
        :param patient: Patient object
        :return: {(prot1, prot2)}
        """
        raise NotImplementedError
        # search strings for protein interaction check if any interacting protein is in patient make sure to delete duplicates. might be sensible to writ Pairs class taking two mutation as input

    def model_pairs_alphafold(self, mutation_a, mutation_b, job_name='', ipath="DB/alphafold_jobs", opath='', padding=100):
        """
        sends two mutations to be modeled in alphafold
        :param mutation_a: Mutation object
        :param mutation_b: Mutation object
        :param name: optional output name
        :param ipath: optional path for fasta files to be sent as jobs
        :param padding: padding on both sides of mutation of modeling sequence
        :param opath: optional output path
        """
        mutation_a.set_ref_seqs_len(padding)
        mutation_b.set_ref_seqs_len(padding)
        seq_1, seq_2 = next(iter((mutation_a.ref_seqs.values()))), next(iter((mutation_b.ref_seqs.values())))
        if (not seq_1) or (not seq_2):
            name = mutation_a.long_name if not seq_1 else mutation_b.long_name
            warnings.warn(f"could not find reference sequence for {name}. Will not predict structure")
            return -1
        name = job_name if job_name else f"{mutation_a.long_name}_{mutation_b.long_name}_{padding}"
        fasta_content = f"{seq_1}:{seq_2}"
        assert len(fasta_content) <= (padding*4)+2, f"{len(fasta_content)} - sequence len appears to be wrong"
        make_fasta(ipath, name, f"{seq_1}:{seq_2}")
        #TODO might want to add option to run the job

    def interactions(self, protein):
        """
        returns list on know interactions
        :param protein: Protein object
        :return: [protein_name]
        """

        return self.protein_nebs[protein.name] if protein.name in self.protein_nebs else []

    '''
    class ProteinScores:
        """
        This class gives a simple API to get all scores related to protein function
        """
        def __init__(self):
            self.distance = 0.0

        @property
        def d(self):
            """
            :return: distance of mutation from catalytic site
            given as mean over all pdbs featuring the mutation
            """
            return self.distance

        @d.setter
        def d(self, value):
            """
            setter for distance score
            :param value: float
            :return: None
            """
            self.distance = value
    '''