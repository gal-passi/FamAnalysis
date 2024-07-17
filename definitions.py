from os.path import dirname, abspath
from os.path import join as pjoin
from hashlib import sha256

hash_url = lambda url: sha256(url).hexdigest()

#  CONTACT INFO CHANGE TO YOUR EMAIL

CONTACT = "gal.passi@mail.huji.ac.il"
HEADERS = {'User-Agent': 'Python {}'.format(CONTACT)}

#  VERBOSE THRESHOLDS

VERBOSE = {'critical': 0, 'program_warning': 1, 'program_progress': 1,
           'thread_warnings': 2, 'thread_progress': 3, 'raw_warnings': 3}
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

#  DIRECTORIES

ROOT_DIR = dirname(abspath(__file__))
DB = 'DB'
DB_PATH = pjoin(ROOT_DIR, DB)
PROTEINS = 'test_p'
PROTEIN_PATH = pjoin(DB_PATH, PROTEINS)
MUTATIONS = 'test_m'
MUTATION_PATH = pjoin(DB_PATH, MUTATIONS)
AFM_PATH = pjoin(DB_PATH, 'AFM')
EVE_PATH = pjoin(DB_PATH, 'EVE')
CPT_PATH = pjoin(DB_PATH, 'CPT')
ESM_PATH = pjoin(DB_PATH, 'ESM')
DB_CHILDREN = ['AFM', 'CPT', 'EVE', 'ESM' 'Family', 'mutations', 'Patients', 'proteins', PROTEINS, MUTATIONS]
CHILDREN_INDEX = {child: i for i, child in enumerate(DB_CHILDREN)}

#  REQUESTS CONSTANTS

TIMEOUT = 10.0
WAIT_TIME = 1.0
RETRIES = 10
RETRY_STATUS_LIST = [429, 500, 502, 503, 504]
DEFAULT_HEADER = "https://"

#  MEMORY CONSTANTS

DEFAULT_RAM_USAGE = 0.65
LOW_MEM_RAM_USAGE = 0.25

#  URLs

UNIPORT_URL = "https://www.uniprot.org/uniprot/"
UNIPORT_QUERY_URL = "https://rest.uniprot.org/uniprotkb/search?"
EBI_PDB_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/"
PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/"
ALPHAFOLD_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v1.pdb"
ALPHAFOLD_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v1.pdb"

# QUERIES

Q_ISOFORMS_PROT = "fields=id&format=tsv&query=gene_exact:{}+AND+organism_id:9606"
Q_ISOFORMS_KEY = "fields=id&format=tsv&query={}+AND+organism_id:9606"
Q_PDBS_UID = "fields=id,xref_pdb&format=tsv&query={}"
#  first protein ref name, second true/false
Q_UID_PROT = "fields=&id&format=tsv&query={}+AND+organism_id:9606+AND+reviewed:{}"
Q_UNIP_ENTERY = "https://rest.uniprot.org/uniprotkb/search?fields=&gene&format=tsv&query={}+AND+organism_id:9606"
Q_UNIP_ENTERY_ALIAS = "https://rest.uniprot.org/uniprotkb/search?fields=&gene&format=tsv&query={}"

# ERRORS & WARNINGS

CON_ERR_EI = "Connection Error in expand_isoforms failed to fetch Uniprot IDs for potein {}"
CON_ERR_FP_1 = "Connection Error in fetch_pdbs for unable to fetch pdbs for Uniprot id: {}\nreturning empty..."
CON_ERR_FP_2 = "Connection Error in fetch_pdbs skipping pdb id {} for Uniprot id {}"
CON_ERR_UFN = "Connection Error in uid_from_name failed to fetch Uniprot IDs for protein {}"
CON_ERR_FUS = "Connection Error in fetch_uniport_sequences while fetching isoforms for {}\nURL: "
CON_ERR_GENERAL = "Connection Error in {} on protein {}"
CON_ERR_DOWNLOAD_SOURCE = "Connection Error while downloading {} data"

DOWNLOAD_INCOMPLETE_WRN = 'File {} is incomplete. Resume download.'
DOWNLOAD_COMPLETE_MSG = 'File {} is complete. Skip download.'
DOWNLOAD_START_MSG = 'File {} does not exist. Start download.'
DOWNLOAD_CORRUPTION_ERR = 'File {} is corrupt. Delete it manually and restart the program.'
DOWNLOAD_NO_HASH_ERR = 'File {} has no hash.'
DOWNLOAD_VALIDATION_MSG = 'File {} is validated.'

#  AMINO ACIDS UTILS

AA_SYN = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS", "I": "ILE",
          "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER",
          "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}
AA_SYN_REV = dict((v, k) for k, v in AA_SYN.items())
AA_TO_INDEX_EVE = {"A": 0, "C": 1, "D": 2, "E": 3, "F": 4, "G": 5, "H": 6, "I": 7, "K": 8, "L": 9, "M": 10,
                   "N": 11, "P": 12, "Q": 13, "R": 14, "S": 15, "T": 16, "V": 17, "W": 18, "Y": 19}
AA_TO_INDEX_ESM = {'K': 0, 'R': 1, 'H': 2, 'E': 3, 'D': 4, 'N': 5, 'Q': 6, 'T': 7, 'S': 8, 'C': 9, 'G': 10,
                   'A': 11, 'V': 12, 'L': 13, 'I': 14, 'M': 15, 'P': 16, 'Y': 17, 'F': 18, 'W': 19}

#  SCORING MODELS

FIRM_SCORE = 'firmScore'  # deprecated
EVE_SCORE = 'eveScore'
EVE_PREDICTION = 'evePrediction'
ESM_SCORE = 'bertScore'
AFM_SCORE = 'afmScore'
EVE_TYPE = 'eveMethod'
AVAILABLE_MODELS = {'EVE', 'ESM', 'AFM'}
MODELS_SCORES = {'EVE': 'eveScore', 'ESM': 'bertScore', 'AFM': 'afmScore', 'FIRM': 'firmScore',
                 'EME_METHOD': 'eveMethod'}
NO_SCORE = -1.0

# ANALYZER CONSTANTS


#  PROTEINS CONSTANTS

PROTEIN_ALIASES = {'LOC100287896': 'LIPT2', 'FPGT-TNNI3K': 'TNNI3K', 'ATPSJ2-PTCD1': 'PTCD1', 'CCL4L1': 'CCL4L2',
                   'PTGDR2': 'CCDC86', '4-SEPT': 'SEPT4'}
NEW_MUTATION_DATA = {'chr': None, 'ref_na': None, 'alt_na': None, 'start': None, 'end': None, AFM_SCORE: NO_SCORE,
                     EVE_SCORE: NO_SCORE, ESM_SCORE: tuple(), EVE_PREDICTION: NO_SCORE, EVE_TYPE: 'no_score'}

#  ALPHA MISSENSE DATA

AFM_PUBLIC_DATA = 'https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz'
AFM_DATA = 'AlphaMissense_aa_substitutions.tsv'
AFM_DATA_PATH = pjoin(AFM_PATH, AFM_DATA)
AFM_DIRECTORY = 'index.json'
AFM_RANGES = 'ranges.json'
AFM_DIRECTORY_PATH = pjoin(AFM_PATH, AFM_DIRECTORY)
AFM_RANGES_PATH = pjoin(AFM_PATH, AFM_RANGES)
AFM_HEADER = 3
AFM_COL_NAMES = ['uniprot_id', 'protein_variant', 'am_pathogenicity', 'am_class']
AFM_ROWSIZE = 32.0
AFM_UID_ROWSIZE = 8.0
AF_ISO_NAME = 'alpha'
AFM_ROWS = 216175351

#  EVE AND CPT DATA

EVE_PUBLIC_DATA = 'https://evemodel.org/api/proteins/bulk/download/?variants=True'
EVE_SINGLE_PROTEIN = 'https://evemodel.org/api/proteins/web_pid/{}/download/?variants=True'
#  to avoid downloading the entire eve dataset insert specific entities e.g BRCA1_HUMAN
EVE_PARTIAL_DOWNLOAD = []
EVE_INDEX_PATH_2 = pjoin(EVE_PATH, 'eve_index.txt')
EVE_EXTENDED_INDEX_PATH = pjoin(EVE_PATH, 'index_unip_full')
EVE_INVERSE_INDEX = pjoin(EVE_PATH, 'eve_reverse_index.txt')
EVE_DATA = 'EVE_all_data'
EVE_DATA_PATH = pjoin(EVE_PATH, EVE_DATA)
EVE_VARIANTS = 'variant_files'
EVE_VARIANTS_PATH = pjoin(EVE_DATA_PATH, EVE_VARIANTS)
EVE_PROT_DOWNLOAD_MSG = "Downloading protein {} from EVE"
EVE_TRUNCATE_TAIL = 10

CPT_DOWNLOAD_BASE = 'https://zenodo.org/records/7954657/files/'
CPT_EVE_DATA = f'{CPT_DOWNLOAD_BASE}CPT1_score_EVE_set.zip?download=1'
CPT_EVE_DATA_HASH = hash_url(CPT_EVE_DATA.encode())
CPT_EVE_DATA_NAME = 'CPT1_score_EVE_set'
CPT_IMPUTE_DATA_1 = f'{CPT_DOWNLOAD_BASE}CPT1_score_no_EVE_set_1.zip?download=1'
CPT_IMPUTE_DATA_1_HASH = hash_url(CPT_IMPUTE_DATA_1.encode())
CPT_IMPUTE_DATA_1_NAME = 'CPT1_score_no_EVE_set_1'
CPT_IMPUTE_DATA_2 = f'{CPT_DOWNLOAD_BASE}CPT1_score_no_EVE_set_2.zip?download=1'
CPT_IMPUTE_DATA_2_HASH = hash_url(CPT_IMPUTE_DATA_2.encode())
CPT_IMPUTE_DATA_2_NAME = 'CPT1_score_no_EVE_set_2'
CPT_IMPUTE_DATA_NAME = 'CPT1_score_no_EVE_set'
CPT_EVE_DATA_PATH = pjoin(CPT_PATH, CPT_EVE_DATA_NAME)
CPT_EVE_IMPUTE_PATH_1 = pjoin(CPT_PATH, CPT_IMPUTE_DATA_1_NAME)
CPT_EVE_IMPUTE_PATH_2 = pjoin(CPT_PATH, CPT_IMPUTE_DATA_2_NAME)
CPT_INGENE_PATH = pjoin(CPT_PATH, CPT_EVE_DATA_NAME)
CPT_EXGENE_PATH = pjoin(CPT_PATH, CPT_IMPUTE_DATA_NAME)
CPT_SCORE_COLUMN = "CPT1_score"
CPT_MUTATION_COLUMN = 'mutant'
CPT_TRUNCATE_TAIL = 7

#  ESM-1b

ESM_VARIANTS_DATA = 'https://huggingface.co/spaces/ntranoslab/esm_variants/resolve/main/ALL_hum_isoforms_ESM1b_LLR.zip'
ESM_DATA = 'ESM_1b_variants'
ESM_DATA_PATH = pjoin(ESM_PATH, ESM_DATA)
ESM_INDEX_PATH = pjoin(ESM_PATH, 'index.json')
ESM_VARIANTS_PATH = pjoin(ESM_DATA_PATH, 'content', 'ALL_hum_isoforms_ESM1b_LLR')
ESM_FILE_SUFFIX = '_LLR.csv'
