from os.path import dirname, abspath
from os.path import join as pjoin

#  CONTACT INFO CHANGE TO YOUR EMAIL

CONTACT = "gal.passi@mail.huji.ac.il"

#  DIRECTORIES

ROOT_DIR = dirname(abspath(__file__))
DB = 'DB'
DB_PATH = pjoin(ROOT_DIR, DB)
PROTEINS = 'test_p'
PROTEIN_PATH = pjoin(DB_PATH, PROTEINS)
MUTATIONS = 'test_m'
MUTATION_PATH = pjoin(DB_PATH, MUTATIONS)

#  URLs

UNIPORT_URL = "https://www.uniprot.org/uniprot/"
UNIPORT_QUERY_URL = "https://rest.uniprot.org/uniprotkb/search?"
EBI_PDB_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/"
PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/"
ALPHAFOLD_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v1.pdb"

#  REQUESTS CONSTANTS

TIMEOUT = 10.0
WAIT_TIME = 1.0
RETRIES = 10
RETRY_STATUS_LIST = [429, 500, 502, 503, 504]
DEFAULT_HEADER = "https://"

#  AMINO ACIDS UTILS

AA_SYN = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS", "I": "ILE",
              "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER",
              "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}
AA_SYN_REV = dict((v, k) for k, v in AA_SYN.items())

#  PROTEINS CHANGE OF NAMES

PROTEIN_ALIASES = {'LOC100287896': 'LIPT2', 'FPGT-TNNI3K': 'TNNI3K', 'ATPSJ2-PTCD1': 'PTCD1', 'CCL4L1': 'CCL4L2',
                'PTGDR2': 'CCDC86', '4-SEPT': 'SEPT4'}

if __name__ == '__main__':
    print(ROOT_DIR)

