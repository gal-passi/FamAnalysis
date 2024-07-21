import re
import os
import pickle
import warnings

import numpy as np
import Analyze
import Protein as P
from definitions import *
from utils import print_if, warn_if



class Mutation:
    """
    This class stores all information about a given mutation
    """
    A_A = "ACDEFGHIKLMNPQRSTVWXY"
    MUTATION_REGEX = rf'p\.(?P<symbol>(?P<orig>[{A_A}]){{1}}(?P<location>[\d]+)(?P<change>[{A_A}]){{1}})'

    def __init__(self, full_desc, protein, dna_data={}, load_only=False, verbose_level=1):
        """
        :param full_desc: string must contain p.{AA}{location}{AA}
        :param protein: protein object to which the mutation belongs
        :param data: dict of mutation data
        :param recover: optional dict containing the recovered data of the object
        """
        self._symbol = None
        self._orig = None
        self._change = None
        self._loc = None
        self._full_desc = None
        self._protein = None
        self._ref_sequences = None  # initialized only upon request
        self._manual_ref = None
        self._pdbs = None  # initialized only upon request
        self._chr, self._start, self._end, self._orig_NA, self._change_NA  = 0, 0, 0, None, None
        self._interface = {}
        self._v = verbose_level
        if not protein or not full_desc:
            raise ValueError("Usage: Mutation(full_desc, protein_name)")

        protein = P.Protein(ref_name = protein, load_only=load_only) if isinstance(protein, str) else protein
        self._directory = os.path.join(MUTATION_PATH, f"{protein.name}_{self.extract_name(full_desc)}.txt")
        if not os.path.exists(self._directory):
            if load_only: raise NameError("Couldn't find mutation in load only mode")
            self._create_new_instance(full_desc, protein, dna_data)
            self._save_obj(self._directory)
        else:
            self._recover_obj(self._directory)

    def __eq__(self, mutatiion):
        """
        equal operator for two mutation opjects
        :param mutatiion: Mutation object
        """
        return (self.name == mutatiion.name) and (self._protein.name == mutatiion._protein.name)

    def __lt__(self, other):
        """
        customize < function will return true if both mutations are from the same protein
        :param other:
        :return:
        """
        if not isinstance(other, Mutation):
            raise NotImplementedError
        return other._protein == self._protein

    def __gt__(self, other):
        """
        customize > function will return true if both mutations are from the same protein
        :param other:
        :return:
        """
        if not isinstance(other, Mutation):
            raise NotImplementedError
        return other._protein == self._protein

    def __hash__(self):
        return hash((self.extended_description, self._protein.name))

    def _create_new_instance(self, full_desc, protein, data):
        """
        sets all atribute of the object
        :param full_desc:
        :param protein:
        :param data:
        :return:
        """
        res = self._descriptor_processor(full_desc)
        self._symbol, self._orig, self._change, self._loc, self._full_desc = \
            res.group('symbol'), res.group('orig'), res.group('change'), int(res.group('location')), full_desc
        if data:
            try:
                self._chr, self._start, self._end, self._orig_NA, self._change_NA = \
                    data['chr'], data['start'], data['end'], data['ref_na'], data['alt_na']
            except KeyError:
               warn_if(self._v, VERBOSE['thread_warnings'], "Incomplete data supplied: data should contain the keys: chr, start, end\n"
                              "the following properties were set to -1: chr, start, end")
               self._chr, self._start, self._end = -1, -1, -1

        self._protein = protein
        self._save_obj(self._directory)

    def _save_obj(self, path):
        data = {'protein': self._protein.name, 'name': self._symbol, 'orig': self._orig,
                'change': self._change, 'location': self._loc, 'full_desc': self._full_desc,
                'chr': self._chr, 'start': self._start, 'end': self._end,
                'orig_NA': self._orig_NA, 'change_NA': self._change_NA, 'manual_ref': self._manual_ref,
                'interface': self._interface}
        with open(path, "wb") as file:
            pickle.dump(data, file)

    def _recover_obj(self, path):
        with open(path, 'rb') as file:
            if os.path.getsize(path) > 0:
                data = pickle.load(file)
                self._symbol, self._orig, self._change, self._loc, self._full_desc,  = \
                    data['name'], data['orig'], data['change'], data['location'], data['full_desc']
                self._interface = data['interface']
                # TODO give optional load for missing data
                self._chr, self._start, self._end, self._orig_NA, self._change_NA  = \
                    data['chr'], data['start'], data['end'], data['orig_NA'], data['change_NA']
                self._protein = P.Protein(ref_name = data['protein'])
                if 'manual_ref' in data.keys():  # TODO this is a temporal fix as not all mut objects have this attribute
                    self._manual_ref = data['manual_ref']
            else:
                raise NameError(f'Data missing or invalid for {path}. \n'
                                f'Delete instance from DB and recreate the object')

    def _descriptor_processor(self, description):
        """
        check the input is compatible with Mutation object constructor
        :param description:
        :return: re.search obj if successful else raise ValueError
        """
        res = re.search(self.MUTATION_REGEX, description)
        if not res:
            raise ValueError("Invalid input valid format of form p.{AA}{location}{AA}")
        return res

    @property
    def _key(self):
        return self.extended_description, self._protein.name

    @property
    def symbol(self):
        return self._symbol

    @property
    def name(self):
        return self._symbol

    @property
    def protein_name(self):
        return self._protein.name

    @property
    def protein(self):
        return self._protein

    @property
    def origAA(self):
        return self._orig

    @property
    def changeAA(self):
        return self._change

    @property
    def loc(self):
        return self._loc

    @property
    def extended_description(self):
        return self.extract_name(self._full_desc)

    @property
    def long_name(self):
        return f"{self._protein.name}_{self._symbol}"

    @property
    def ref_seqs(self):
        """
        returns a dictionary of 10 AA res_seq in which the mutation is found.
        the sequences are obtained from all known isoforms of the Protein obj of which the mutation is found.
        :return: {isoform: ref_seq}
        """
        if not self._ref_sequences:
            self._find_reference_sequences()
        if self._manual_ref:
            self._ref_sequences['manual_ref'] = self._manual_ref
        return self._ref_sequences

    def set_ref_seqs_len(self, n):
        """
        updates reference sequences property to the desired length
        :param n: the new references will be of length 2n
        :return:
        """
        self._find_reference_sequences(padding=n)

    def set_manual_ref_sequence(self, seq):
        """
        updates reference sequences property to the desired length
        :param n: the new references will be of length 2n
        :return:
        """
        self._ref_sequences = {'manual':  seq}


    @property
    def pdbs(self):
        """
        :return: {isoform: pdbs_id} or {} if none are found
        """
        if not self._pdbs:
            self._find_pdbs()
        return self._pdbs

    @property
    def raw_pdbs(self):
        """
        :return: {pdb_id} or {} if none are found
        """
        pdbs = set()
        for tuples in self.pdbs.values():
            if tuples:
                for id, _ in tuples:
                    if id == 'manual':
                        continue
                    pdbs.add(id.lower())
        return pdbs

    @property
    def interface(self):
        if self._interface: return self._interface
        else: return -1
        #TODO fix should not calculate by default
        #if not self.raw_pdbs: return {}
        #self._interface = self.calc_interface()
        #return self._interface

    @property
    def chr(self):
        return self._chr

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    @property
    def ref_na(self):
        return self._orig_NA

    @property
    def alt_na(self):
        return self._change_NA

    @property
    def isoforms(self):
        """
        returns a set of all valid isoforms
        """
        isoforms = set(self.ref_seqs.keys())
        isoforms.discard('manual_ref')
        if isoforms:
            return isoforms
        elif self._manual_ref:
            isoforms = [iso for iso, seq in self._protein.isoforms.items() if self._manual_ref in seq]
            return set(isoforms)
        return set()

    @property
    def mentioned(self):
        """
        returns whether the protein was mentioned in litriture
        """
        return self._protein.mentioned

    @property
    def for_bert_location(self, offset=1):
        return self.origAA + str(self.loc-offset) + self.changeAA

    @property
    def esm_score(self):
        """
        :return: esm-2v score
        """
        try:
            return self._protein.muts[self.extended_description][MODELS_SCORES['ESM']]
        except KeyError:
            warn_if(self._v, VERBOSE['thread_warnings'], f"eveScore not initialized for mutation - try to add mutation again")
            warn_if(self._v, VERBOSE['raw_warnings'], f"\n{self._protein.muts}")
            return None


    @property
    def has_esm(self):
        return self.esm_score is not None

    @property
    def eve_score(self):
        try:
            return self._protein.muts[self.extended_description][MODELS_SCORES['EVE']]
        except KeyError:
            warn_if(self._v, VERBOSE['thread_warnings'], f"eveScore not initialized for mutation - try to add mutation again")
            warn_if(self._v, VERBOSE['raw_warnings'], f"\n{self._protein.muts}")
            return -1

    @property
    def has_eve(self):
        return self.eve_score not in [None, -1]

    @property
    def afm_score(self):
        return self._protein.muts[self.extended_description][AFM_SCORE]

    @property
    def has_afm(self):
        return (self.afm_score is not None) and (self.afm_score != NO_SCORE)

    @property
    def eve_prediction(self):
        return self._protein.muts[self.extended_description][EVE_PREDICTION]

    @property
    def eve_type(self):
        return self._protein.muts[self.extended_description][EVE_TYPE]

    @property
    def firm_score(self):
        s = self._protein.muts[self.extended_description]['firmScore']
        if isinstance(s, float) or isinstance(s, int):
            return s
        return -1

    @property
    def consensus_score(self):
        return self._protein.muts[self.extended_description]['consensusScore']

    @property
    def rank_score(self):
        return self._protein.muts[self.extended_description]['rank']

    @property
    def has_firm(self):
        return self.firm_score not in {-1, 'asr', 'err'}

    @property
    def alph_conf(self):
        return self.protein.mutations[self.extended_description]['AlphaFoldConfidence']

    @property
    def n_scores(self):
        return float(self.has_firm + self.has_eve + self.has_esm)


    @staticmethod
    def extract_name(description):
        """
        a static method to extract mutation reference name out of of free text
        :param description: str general description
        :return: string mutation name
        """
        aa = "ACDEFGHIKLMNPQRSTVWXY"
        mutation_regex = rf'(?P<symbol>p\.(?P<orig>[{aa}]){{1}}(?P<location>[\d]+)(?P<change>[{aa}]){{1}})'
        res = re.search(mutation_regex, description)
        if res:
            return res.group('symbol')
        else:
            raise ValueError

    def _find_reference_sequences(self, seq="", padding=5):
        """
        :param protein: Protein object
        :param seq: optional specific sequence to search in
        :param padding: optional sets the length of the reference sequence default is 10
        :return: None sets reference obj
        """
        search_loc = self._loc - 1
        references = {}
        found = set()
        start = 0 if (search_loc - padding < 0) else (search_loc - padding)  # avoids out of range
        if (seq != ''):
            try:
                if seq[search_loc] == self._orig:
                    ref = seq[start: search_loc + padding]
                    found.add(ref)
                    references["Unknown_iso"] = ref

            except IndexError:
                pass

        for iso, seq in self._protein.isoforms.items():
            try:
                if seq[search_loc] == self._orig:
                    ref = seq[start: search_loc + padding]
                    if ref not in found:
                        found.add(ref)
                        references[iso] = ref
            except IndexError:
                pass

        self._ref_sequences = references

    def _find_pdbs(self):
        """
        finds all pdbs relevent to the a_a in location
        ":param location: int location to search
        :param a_a: a_a in location
        :param padding: padding for reference seq optional
        :return: dict {isoform: pdb_ids}
        """
        self._pdbs = {iso: self._find_relevent_pdbs(ref) for iso, ref in self.ref_seqs.items()}

    def _find_relevent_pdbs(self, reference_sequence):
        """
        searches know pdbs that contain the reference sequence
        sequences too short will result in non-specific results
        :param reference_sequence: AA sequence
        :return: list of relevent pdb-ids
        """
        return [(p_id, seq.find(reference_sequence) + 5) for p_id, seq in self._protein.pdbs.items() if (reference_sequence in seq) and (p_id != 'manual')]

    def has_pdbs(self):
        for key, value in self.pdbs.items():
            if value:
                return True
        return False

    def calc_interface(self):
        return Analyze.ProteinAnalyzer().score_distances(self)

    def min_interface(self):
        if not self.interface: return np.inf
        if isinstance(self.interface, dict):
            return np.amin(np.array(list(self.interface.values())).astype(np.float))
        else:
            return self.interface

    def print_pdbs(self):
        s = ''
        for pdbs in self.pdbs.values():
            for pdb, _ in pdbs:
                s += pdb + " "
        return s

    def download_pdbs(self, path='DB/pdbs'):
        """
        downloads all known pdbs of the mutation
        :param path: path to save pdb in
        :return:
        """
        for pdbs in self.pdbs.values():
            for pdb, _ in pdbs:
                self._unip.download_pdb(pdb, path)

    def update_score(self, model, score, eve_type='', esm_type=''):
        """
        updates mutation model scores
        :param model: str one of: EVE | ESM | AFM
        :param score: float score to update
        :param esm_type: specification of model used for eveModel score
        :param eve_type: specification of inference used in esm
        :return:
        """
        assert model in AVAILABLE_MODELS, f"model must be one of {' '.join(AVAILABLE_MODELS)}"
        assert isinstance(score, (float, int, type(None))), 'score must be either float or None'
        prot = self.protein.reload()
        prot.muts[self.extended_description][MODELS_SCORES[model]] = score
        if model == 'EVE' and eve_type:
            prot.muts[self.extended_description][MODELS_SCORES['EVE_METHOD']] = eve_type
        if model == 'ESM' and esm_type:
            prot.muts[self.extended_description][MODELS_SCORES['ESM_METHOD']] = esm_type
        prot._update_DB(pjoin(prot.directory, prot.MUTS), prot.muts, mode='pickle')

    def esm_scores_from_inference(self, how='all'):
        """
        This method was deprecated - used for direct inference from esm
        :param how: one of ['all', 'sum', 'avg', 'min', 'max']
        :return:
        """
        warnings.warn('This method is no longer supported use esm_score instead')
        return None
        scores = self._protein.muts[self.extended_description]['bertScore']
        if not scores or type(scores) != tuple:
            return 0
        if how == 'all':
            return scores
        if how == 'sum':
            return sum(scores)
        if how == 'avg':
            return sum(scores)/float(len(scores))
        if how == 'min':
            return min(scores)
        if how == 'max':
            return max(scores)

    def print_status(self):
        print(f"protein name:           {self.protein_name}")
        print(f"mutation symbol:        {self.name}")
        print(f"esm score (min):       {round(self.esm_score, 3)}")
        try:
            print(f"firm score:             {round(self.firm_score, 3)}")
        except TypeError:
            print(f"firm score:             {self.firm_score}")
        print(f"eve score & status:     {round(self.eve_score, 3)}, {self.eve_prediction}")
        print(f"interface score:        {self.min_interface()}")
        print(f"consensus_score:        {self.consensus_score}")
        print(f"rank:                   {self.rank_score}")
        print(f"litereture mention:     {self.mentioned}")
        print(f"matching pdbs:          {self.print_pdbs()}")
        print(f"proteins pdbs:          {list(self.protein.pdbs.keys())}")
        print(f"alphafold conf:         {self.alph_conf}")
        print(f"known interactions:     {self.protein.interactions}")