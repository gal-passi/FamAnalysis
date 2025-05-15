import os.path
import json
import pickle
import re
from datetime import datetime
import Mutation
from Connections import Uniport
from Connections import EntrezApi as entrez
from definitions import *
from utils import print_if, warn_if
from copy import deepcopy


class Protein:
    """
    This class creates and retrieves local proteins DBs
    """
    MAIN_SEQ = "sequence.txt"
    ISOFORMS = "isoforms.txt"
    UIDS = "uids.txt"
    PDBS = "pdbs.txt"
    MUTS = "Mutations.txt"
    REF_SEQ = "reference.txt"
    BACKUP_PATH = "backup.txt"

    def __init__(self, ref_name='', uniport_id="", ncbi="", load_only=False, verbose_level=1):
        """
        Constructor for Protein Class. Must either contain free_text or ref_name
        :param ref_name: optional initializor using protein ref name i.e. LRRK2
        :param uniport_id: str optional uniport id / primary accession. if not given the id will be extracted using the protein ref name
        :param ncbi: if given will augment protein isoforms with the isoforms matching the id str | set | list
        :param add_mut: bool if set to true the mutation in the free text will be added to the protein's mutations dict
        """
        if (not ref_name) and (not uniport_id):
            raise ValueError("ERROR: Class initialization must contain either protein reference name "
                             "or uniprot id \n" "USAGE: Protein('CUL7') | Protein(uniprot_id = Q14999)")
        self._v = verbose_level
        self._unip = Uniport(verbose_level=verbose_level)
        self._ncbi_ids = set()
        if isinstance(ncbi, set):
            self._ncbi_ids |= ncbi
        elif isinstance(ncbi, list):
            self._ncbi_ids |= set(ncbi)
        elif isinstance(ncbi, str) and ncbi:
            self._ncbi_ids.add(ncbi)
        ref_name = PROTEIN_ALIASES[ref_name] if ref_name in PROTEIN_ALIASES else ref_name
        self._name = uniport_id if not ref_name else ref_name
        self.directory = os.path.join(PROTEIN_PATH, self._name)
        if not os.path.exists(self.directory):
            if not load_only:
                print_if(self._v, VERBOSE['thread_progress'], f"Creating new protein {ref_name}:")
                self.create_new_entity(uniport_id)
            else:
                raise NameError("Couldn't find protein in load only mode")

        else:
            print_if(self._v, VERBOSE['thread_progress'], f"Protein {ref_name} was loaded")
            self._Uids, self.isoforms, self.pdbs, self.muts = self.load_protein()

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        if not isinstance(other, Protein):
            raise NotImplementedError
        return self.name == other.name

    def reload(self):
        """
        reload object from memory
        :return:
        """
        self._Uids, self.isoforms, self.pdbs, self.muts = self.load_protein()
        return self

    @property
    def name(self):
        """
        return proteins symbol
        :return: str symbol
        """
        return self._name

    @property
    def mutations(self):
        return list(self.muts.keys())

    def generate_mutations(self):
        """
        generator of all self Mutation objects
        :return:
        """
        for mut_desc in self.mutations:
            yield Mutation.Mutation(mut_desc, self)

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
    def aliases(self):
        ret = set(self._Uids['aliases'])
        ret = ret - {''}  # remove empty strings
        return ret

    def entry_name(self, all=False):
        """
        :return: list of uniprot accessions
        """
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
        :param mut: string p.{AA}{location}{AA}
        :return: {isoform: [pdb_pathes]}
        """
        if mut not in self.muts:
            print_if(self._v, VERBOSE['thread_warnings'],
                     f"{mut} not found make sure {mut} is added to {self.name} using self.add_mut")
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
        uids = self._unip.uid_from_name(self.name)
        if uniprot_id and (uniprot_id not in uids['reviewed']):
            uids['reviewed'].append(uniprot_id)
        self._Uids = uids
        if self.Uid == '' and not self._ncbi_ids:
            raise NameError(f'Unable to find uniport id for protein {self._name}\n'
                            f'Protein will not be created')
        try:
            os.mkdir(self.directory)
        except OSError:
            raise Exception(f"ERROR: FAILED TO CREATE FILE {self.directory}")

        with open(os.path.join(self.directory, self.UIDS), "w") as file:
            file.write(json.dumps(uids))

        isoforms = self._unip.fetch_uniport_sequences(self.Uid)
        alphafold_seq = self._unip.alpha_seq(self)
        ncbi = {} if not self.ncbi_ids else entrez().fetch_NCBI_sequences(self.ncbi_ids)
        isoforms = {**isoforms, **alphafold_seq, **ncbi}
        self.isoforms = isoforms
        with open(os.path.join(self.directory, self.ISOFORMS), "w") as file:
            file.write(json.dumps(isoforms))
        pdbs = self._unip.fetch_pdbs(prot=self)
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

    def add_isoform(self, id, dataset):
        """
        :param id: str
        :param dataset: str either ncbi | uniprot
        :return:
        """
        if not isinstance(id, str) or not dataset in ['ncbi', 'uniprot']:
            raise ValueError
        if dataset == 'ncbi':
            self.isoforms |= entrez().fetch_NCBI_sequences(id)
        elif dataset == 'uniprot':
            self.isoforms |= self._unip.fetch_uniport_sequences(id)
        self._update_DB(os.path.join(self.directory, self.ISOFORMS), self.isoforms, mode='json')
        return self.reload()

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
            warn_if(self._v, VERBOSE['thread_warnings'], f"Unable to find mutation name in reference:\n{description}")
            return

        if name in self.muts:
            return

        try:
            Mutation.Mutation(description, self, dna_data, verbose_level=self._v)  # check if given mutation is valid
        except ValueError:
            warn_if(self._v, VERBOSE['thread_warnings'], f"Unable to find mutation in reference:\n{description}")
            return

        data = deepcopy(NEW_MUTATION_DATA)
        if dna_data:
            try:
                data['chr'] = dna_data['chr']
                data['ref_na'] = dna_data['ref_na']
                data['alt_na'] = dna_data['alt_na']
                data['start'] = dna_data['start']
                data['end'] = dna_data['end']
            except KeyError:
                warn_if(self._v, VERBOSE['thread_warnings'],
                        f"DNA data supplied for {description} in wrong format skipping:\n{dna_data}")
                data = deepcopy(NEW_MUTATION_DATA)

        self.muts[name] = data
        self._update_DB(os.path.join(self.directory, self.MUTS), self.muts, mode='pickle')
        print_if(self._v, VERBOSE['thread_progress'], f"added {description} to {self._name}")

    def find_relevent_pdbs(self, reference_sequence):
        """
        searches know pdbs that contain the reference sequence
        sequences too short will result in non-specific results
        :param ref: AA sequence
        :return: list of relevent pdb-ids
        """
        return [p_id for p_id, seq in self.pdbs.items() if reference_sequence in seq]
