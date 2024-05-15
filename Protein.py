import os.path
import json
import pickle
import re
from datetime import datetime
import warnings
import Mutation
import Analyze
from Connections import Uniport

PATH = r"sequence.txt"
PAD = 6


class Protein:
    """
    This class creates and retrieves local proteins DBs
    """
    PROTEIN_DB = "DB/proteins"
    MAIN_SEQ = "sequence.txt"
    ISOFORMS = "isoforms.txt"
    UIDS = "uids.txt"
    PDBS = "pdbs.txt"
    MUTS = "Mutations.txt"
    REF_SEQ = "reference.txt"
    BACKUP_PATH = "backup.txt"
    UNIPORT_URL = "http://www.uniprot.org/uniprot/"
    UNIQUE_NAMES = {'LOC100287896': 'LIPT2', 'FPGT-TNNI3K': 'TNNI3K', 'ATPSJ2-PTCD1': 'PTCD1', 'CCL4L1': 'CCL4L2', 'PTGDR2': 'CCDC86'}

    def __init__(self, free_text='', ref_name='', uniport_id="", ncbi="", load_only=False):
        """
        Constructor for Protein Class. Must either contain free_text or ref_name
        :param free_text: optional initialization using free text i.e. LRRK2:NM_198578.1:exon48:c.G7134C:p.K2378N
        :param ref_name: optional initializor using protein ref name i.e. LRRK2
        :param uniport_id: int optopnal uniport id. if not given the id will be extracted using the protein ref name
        :param ncbi: int optional ncbi id if not given the id will be extracted using the protein ref name
        :param add_mut: bool if set to true the mutation in the free text will be added to the protein's mutations dict
        """
        if (free_text == '') and (ref_name == ''):
            raise ValueError("ERROR: Class initialization must contain either protein reference name "
                             "or free description text\n" "USAGE: Protein('CUL7') | "
                             "Protein('CUL7:NM_001168370:exon25:c.G4970A:p.R1657Q')")
        self._unip = Uniport()
        if free_text != '':
            data = re.search(rf'([A-Z\d]*):(NM_[\d]*)[\.]?[\d]*:', free_text)
            ref_name = data.group(1) if data is not None else ref_name
            if ncbi != "":
                ncbi = ncbi
            elif data is not None:
                ncbi = data.group(2)
            else:
                ncbi = ""

        if ref_name == '' and not uniport_id:
            raise ValueError("ERROR: Unable to extract protein reference name, check input\n"
                             "USAGE: Protein('CUL7') | Protein('CUL7:NM_001168370:exon25:c.G4970A:p.R1657Q')")
        self.ncbi_id = ncbi  #TODO not implemented well not sure if even needed
        ref_name = self.UNIQUE_NAMES[ref_name] if ref_name in self.UNIQUE_NAMES else ref_name
        self._name, self.directory = ref_name, os.path.join(self.PROTEIN_DB, ref_name)
        if not os.path.exists(self.directory):
            if not load_only:
                print(f"Creating new protein {ref_name}:")
                self.create_new_entity(uniport_id)
            else:
                raise NameError("Couldn't find protein in load only mode")

        else:
            self._Uids, self.isoforms, self.pdbs, self.muts = self.load_protein()

        #TODO this data should be saved and loaded not re-calculated
        #self.known_oncogene, self.mentioned_in = Analyze.ProteinAnalyzer.search_litriture(self.name)
        #self.evemodel = Analyze.ProteinAnalyzer.search_eve_record(self.Uid)

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        if not isinstance(other, Protein):
            raise NotImplementedError
        return self.name == other.name


    @property
    def name(self):
        """
        return proteins symbol
        :return: str symbol
        """
        return self._name

    @property
    def mutations(self):
        return self.muts

    @property
    def Uid(self):
        """
        :return: First uid found in records - will prefer from reviewed source
        """
        if len(self._Uids['reviewed']) > 0:
            return self._Uids['reviewed'][0]
        elif len(self._Uids['non_reviewed']) > 0:
            return self._Uids['non_reviewed'][0]
        return ''


    @property
    def mentioned(self):
        """
        returns whethere the protein was mentioned in litriture
        """
        return Analyze.ProteinAnalyzer()._search_litriture(self, bool=True)

    @property
    def interactions(self):
        return Analyze.ProteinAnalyzer().interactions(self)

    def entery_name(self, all=False):
        return self._Uids['main_entery'] if not all else self._Uids['all_enteries']

    def all_uids(self):
        return self._Uids

    def has_mutation(self, mutation):
        """
        Mutation obj
        :param mutation:
        :return:
        """
        return mutation.name in self.muts

    def pdb_paths(self, mut):
        """
        return a dict of all pdbs relevent to the mut
        :param mut: strin p.{AA}{location}{AA}
        :return: {isoform: [pdb_pathes]}
        """
        if mut not in self.muts:
            print("Mutation not found please add mutation using the add_mut function and rerun")
            return []
        return {iso: [os.path.join(self.directory, mut, pdb + ".ent") for pdb in self.muts[mut]["pdbs"][iso]]
                for iso in self.muts[mut]["pdbs"]}

    def load_protein(self):
        """
        loads protein from DB
        :param ref_name: refrence name of protein as saved in DB
        :return:
        """
        return self._read_uids(), self._read_isoforms(), self._read_pdbs(), self._read_mutations()
        #return self._read_isoforms(), self._read_pdbs(), self._read_mutations()

    def _read_uids(self):
        return self.read_file(os.path.join(self.directory, self.UIDS))

    def _read_isoforms(self):
        return self.read_file(os.path.join(self.directory, self.ISOFORMS))

    def _read_pdbs(self):
        return self.read_file(os.path.join(self.directory, self.PDBS))

    def _read_mutations(self):
        """
        restores Protein's mutation dictionary.
        the data is saved as a set of mutation names restores a dict {mutation_name: Mutation obj}
        :return:
        """
        return self.read_file(os.path.join(self.directory, self.MUTS), mode=3)

    @staticmethod
    def read_file(path, mode=1):
        """
        Reads amino acids sequence - returns a string with no whitespaces
        :param path:
        :param mode: file mode default is 1 - json format, 2 - plain text, 3 - pickle format
        :return: the read file is mode = 1 returned as dict else as text
        """
        read_mode = "rb" if mode == 3 else "r"
        with open(path, read_mode) as file:
            if mode == 1:
                return json.loads(file.read())
            elif mode == 2:
                return json.load(file)
            elif mode == 3:
                if os.path.getsize(path) > 0:
                    return pickle.load(file)
                else:
                    return set()
            else:
                raise ValueError("Invalid mode should either 1 or 2")

    def create_new_entity(self, uniprot_id=""):
        """
        creats a new protein entity
        :param name:
        :param uniprot_id: optional if given will look for protein data using the given id
        :return:
        """
        uids = {'reviewed': self._unip.uid_from_name(self.name, reviewed=True, all=True),
                'non_reviewed': self._unip.uid_from_name(self.name, reviewed=False, all=True),
                'main_entery': self._unip.entery_name(self),
                'all_enteries': self._unip.entery_name(self, all_results=True)}
        if uniprot_id:
            uids['reviewed'].append(uniprot_id)
        self._Uids = uids
        if self.Uid == '' and not self.ncbi_id:
            raise NameError(f'Unable to find uniport id for protein {self._name}\n'
                            f'Protein will not be created')
        try:
            os.mkdir(self.directory)
        except OSError:
            raise Exception(f"ERROR: FAILED TO CREATE FILE {self.directory}")

        with open(os.path.join(self.directory, self.UIDS), "w") as file:
            file.write(json.dumps(uids))

        isoforms = self._unip.fetch_uniport_sequences(self.Uid) #TODO might want to expend isoform selection
        ncbi = {} if not self.ncbi_id else self._unip.fatch_all_NCBIs(self.ncbi_id)
        isoforms = {**isoforms, **ncbi}
        self.isoforms = isoforms
        with open(os.path.join(self.directory, self.ISOFORMS), "w") as file:
            file.write(json.dumps(isoforms))

        pdbs = self._unip.fetch_pdbs(prot=self, verbal=True)
        self.pdbs = pdbs
        with open(os.path.join(self.directory, self.PDBS), "w") as file:
            file.write(json.dumps(pdbs))

        self.muts = {}
        with open(os.path.join(self.directory, self.MUTS), "wb") as file:
            pickle.dump(dict(), file)

        open(os.path.join(self.directory, self.BACKUP_PATH), "w").close()

    def _update_DB(self, to_update, data, mode='json'):
        """
        updates textual DB entity
        :param to_update: entity to update
        :param data: variable to insert to DB
        :param mode: either json (default) or pickle for the mode in which the variable needs to be saved
        :return: bool
        """
        # TODO backup might not work on binary files used in pickle mode
        open_mode = "r+" if mode == 'json' else 'rb+'
        with open(to_update, open_mode) as file:
            backup = file.read()
            file.truncate(0)

        open_mode = 'w' if mode == 'json' else 'wb'
        with open(to_update, open_mode) as file:
            if mode == 'json':
                file.write(json.dumps(data))
            elif mode == 'pickle':
                pickle.dump(data, file)

        with open(os.path.join(self.directory, self.BACKUP_PATH), "a") as file:
            str_time = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
            file.write(f"changed {to_update} on {str_time}. Erased: \n {backup} \n")

    def add_mut(self, description, dna_data={}):
        """
        adds mutation to protein DB
        :param description: mutation discription must contain p.{AA}{location}{AA}
        :param dna_data: optinal dict of mutation dna_data in format: {'chr': , 'start': , 'end' : , 'ref_na': , 'alt_na': }
        :return: dict {"details" : (full description, AA, int(location), AA, [[ref_seq, iso]]),
                                          "pdbs" : {iso: [pdbs]}
        """
        try:
            name = Mutation.Mutation.extract_name(description)
        except ValueError:
            warnings.warn(f"Unable to find mutation in reference:\n{description}")
            return

        if name in self.muts:
            return

        try:
            Mutation.Mutation(description, self, dna_data)  # check if given mutation is valid
        except ValueError:
            warnings.warn(f"Unable to find mutation in reference:\n{description}")
            return

        data = {'chr': None, 'ref_na': None, 'alt_na': None, 'start': None, 'end': None, 'firmScore': -1.0,
                'eveScore': -1.0, 'bertScore': tuple(), 'evePrediction': -1.0}
        if dna_data:
            try:
                data['chr'] = dna_data['chr']
                data['ref_na'] = dna_data['ref_na']
                data['alt_na'] = dna_data['alt_na']
                data['start'] = dna_data['start']
                data['end'] = dna_data['end']
            except KeyError:
                warnings.warn(f"DNA data supplied in wrong format skipping:\n{dna_data}")
                data = {'chr': None, 'ref_na': None, 'alt_na': None, 'start': None, 'end': None, 'firmScore': -1.0,
                        'eveScore': -1.0, 'bertScore': tuple(), 'evePrediction': -1.0}

        self.muts[name] = data
        self._update_DB(os.path.join(self.directory, self.MUTS), self.muts, mode='pickle')
        print(f"added {description} to {self._name}")

    def find_relevent_pdbs(self, reference_sequence):
        """
        searches know pdbs that contain the reference sequence
        sequences too short will result in non-specific results
        :param ref: AA sequence
        :return: list of relevent pdb-ids
        """
        return [p_id for p_id, seq in self.pdbs.items() if reference_sequence in seq]



#return {mutation: Mutation.Mutation(mutation, self) for mutation in mutations}