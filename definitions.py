from os.path import dirname, abspath
from os.path import join as pjoin
from hashlib import sha256
import torch
import re
from pathlib import Path

hash_url = lambda url: sha256(url).hexdigest()

#  CONTACT INFO CHANGE TO YOUR EMAIL

CONTACT = ""
HEADERS = {'User-Agent': 'Python {}'.format(CONTACT)}
HUGGINGFACE_TOKEN = "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"  # obtain read permission token

#  DEVICE

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'

#  VERBOSE THRESHOLDS

VERBOSE = {'critical': 0, 'program_warning': 1, 'program_progress': 1,
           'thread_warnings': 2, 'thread_progress': 3, 'raw_warnings': 3}
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

#  DIRECTORIES

ROOT_DIR = dirname(abspath(__file__))
DB = 'DB'
DB_PATH = pjoin(ROOT_DIR, DB)
PROTEINS = 'proteins'
PROTEIN_PATH = pjoin(DB_PATH, PROTEINS)
MUTATIONS = 'mutations'
MUTATION_PATH = pjoin(DB_PATH, MUTATIONS)
AFM_PATH = pjoin(DB_PATH, 'AFM')
EVE_PATH = pjoin(DB_PATH, 'EVE')
EVE_PARQUET_DIR = Path(EVE_PATH) / "parquet_files"
CPT_PATH = pjoin(DB_PATH, 'CPT')
ESM_PATH = pjoin(DB_PATH, 'ESM')
DB_CHILDREN = ['AFM', 'CPT', 'EVE', 'ESM', 'Family', 'Patients', PROTEINS, MUTATIONS]
CHILDREN_INDEX = {child: i for i, child in enumerate(DB_CHILDREN)}
EVE_SETUP_INDEX_DIR = 'eve_index'

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
ENTREZ_API_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
ENTREZ_AA_RET = 'fasta_cds_aa'
ENTREZ_NA_RET = 'fasta_cds_na'
ENTREZ_SEARCH_URL = ENTREZ_API_URL + "esearch.fcgi?db=nuccore&term={}&usehistory=y"
ENTREZ_SEQ_FROM_KEYS = ENTREZ_API_URL + \
                       "efetch.fcgi?db=nuccore&query_key={}&WebEnv={}&rettype={}&retmode=text"

ENTREZ_WEBENV_RE = re.compile(r'\<WebEnv\>(\S+)\<\/WebEnv\>')
ENTREZ_QUERYKEY_RE = re.compile(r'\<QueryKey\>(\d+)\<\/QueryKey\>')

# QUERIES

Q_ISOFORMS_PROT = "fields=id&format=tsv&query=gene_exact:{}+AND+organism_id:9606"
Q_ISOFORMS_KEY = "fields=id&format=tsv&query={}+AND+organism_id:9606"
#  first protein ref name, second true/false
Q_UID_PROT = "fields=&id&format=tsv&query={}+AND+organism_id:9606+AND+reviewed:{}"
Q_UID_PROT_ALL = "fields=&gene&format=tsv&query={}+AND+organism_id:9606"
Q_UNIP_ENTERY = "https://rest.uniprot.org/uniprotkb/search?fields=&gene&format=tsv&query={}+AND+organism_id:9606"
Q_UNIP_ENTERY_ALIAS = "https://rest.uniprot.org/uniprotkb/search?fields=&gene&format=tsv&query={}"
Q_UNIP_ALL_ISOFORMS = UNIPORT_QUERY_URL + "&format=fasta&query=" \
                                          "(accession:{}+AND+is_isoform:true)+OR+(accession:{}+AND+is_isoform:false)"
Q_UNIP_ALL_PDBS = UNIPORT_QUERY_URL + "fields=id,xref_pdb&format=tsv&query={}"

UIDS_COL_IDX = 0
REVIEWED_COL_IDX = 2
GENE_NAME_COL_IDX = 4
CHUNK_SIZE = 10
UNIP_REVIEWED = 'reviewed'

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

A_A = "ACDEFGHIKLMNPQRSTVWXY"
N_AA = 20
AA_SYN = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS", "I": "ILE",
          "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER",
          "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}
AA_SYN_REV = dict((v, k) for k, v in AA_SYN.items())
AA_TO_INDEX_EVE = {"A": 0, "C": 1, "D": 2, "E": 3, "F": 4, "G": 5, "H": 6, "I": 7, "K": 8, "L": 9, "M": 10,
                   "N": 11, "P": 12, "Q": 13, "R": 14, "S": 15, "T": 16, "V": 17, "W": 18, "Y": 19}
AA_TO_INDEX_ESM = {'K': 0, 'R': 1, 'H': 2, 'E': 3, 'D': 4, 'N': 5, 'Q': 6, 'T': 7, 'S': 8, 'C': 9, 'G': 10,
                   'A': 11, 'V': 12, 'L': 13, 'I': 14, 'M': 15, 'P': 16, 'Y': 17, 'F': 18, 'W': 19}

STOP_CODONS = ['tag', 'taa', 'tga']
STOP_AA = '_'
CODON_LENGTH = 3
CODON_TRANSLATOR = {'ata': 'I', 'atc': 'I', 'att': 'I', 'atg': 'M', 'aca': 'T',
                    'acc': 'T', 'acg': 'T', 'act': 'T', 'aac': 'N', 'aat': 'N',
                    'aaa': 'K', 'aag': 'K', 'agc': 'S', 'agt': 'S', 'aga': 'R',
                    'agg': 'R', 'cta': 'L', 'ctc': ' L', 'ctg': 'L', 'ctt': 'L',
                    'cca': 'P', 'ccc': 'P', 'ccg': 'P', 'cct': 'P', 'cac': 'H',
                    'cat': 'H', 'caa': 'Q', 'cag': 'Q', 'cga': 'R', 'cgc': 'R',
                    'cgg': 'R', 'cgt': 'R', 'gta': 'V', 'gtc': 'V', 'gtg': 'V',
                    'gtt': 'V', 'gca': 'A', 'gcc': 'A', 'gcg': 'A', 'gct': 'A',
                    'gac': 'D', 'gat': 'D', 'gaa': 'E', 'gag': 'E', 'gga': 'G',
                    'ggc': 'G', 'ggg': 'G', 'ggt': 'G', 'tca': 'S', 'tcc': 'S',
                    'tcg': 'S', 'tct': 'S', 'ttc': 'F', 'ttt': 'F', 'tta': 'L',
                    'ttg': 'L', 'tac': 'Y', 'tat': 'Y', 'taa': '_', 'tag': '_',
                    'tgc': 'C', 'tgt': 'C', 'tga': '_', 'tgg': 'W'}

#  SCORING MODELS

FIRM_SCORE = 'firmScore'  # deprecated
EVE_SCORE = 'eveScore'
EVE_PREDICTION = 'evePrediction'
ESM_SCORE = 'bertScore'
ESM3_SCORE = 'esm3Score'
ESM_TYPE = 'esmMethod'
ESM3_TYPE = 'esm3Method'
AFM_SCORE = 'afmScore'
EVE_TYPE = 'eveMethod'
DS_RANK = 'DSRank'
AVAILABLE_MODELS = {'EVE', 'ESM', 'AFM', 'ESM3', 'ESM1B', 'ESM1b'}
AVAILABLE_SCORES = AVAILABLE_MODELS | {'DS'}
MODELS_SCORES = {'EVE': EVE_SCORE, 'ESM': ESM_SCORE, 'ESM1B': ESM_SCORE, 'ESM1b': ESM_SCORE, 'ESM3': ESM3_SCORE,
                 'AFM': AFM_SCORE, 'FIRM': FIRM_SCORE, 'EVE_METHOD': EVE_TYPE, 'ESM_METHOD': ESM_TYPE, 'DS': DS_RANK,
                 'ESM3_METHOD': ESM3_TYPE, 'ESM1b_METHOD': ESM_TYPE, 'ESM_METHOD1B': ESM_TYPE}
NO_SCORE = -1.0
NO_TYPE = 'no_score'

#  PROTEINS CONSTANTS

PROTEIN_ALIASES = {'LOC100287896': 'LIPT2', 'FPGT-TNNI3K': 'TNNI3K', 'ATPSJ2-PTCD1': 'PTCD1', 'CCL4L1': 'CCL4L2',
                   'PTGDR2': 'CCDC86', '4-SEPT': 'SEPT4', '4-Sep': 'SEPT4','TPTEP2-CSNK1E': 'CSNK1E', 'LOC101928841': 'ADPRHL1',
                   'LOC102724159': 'PWP2', 'LOC100421372': 'MSANTD7', 'LOC102724488': 'SYT15B', 'LOC645177': 'IRAG2',
                   'CBSL': 'CBS_HUMAN'}
REMOVED_PROTEINS = ['LOC100996842', 'GLRA4', 'LOC390877', 'NM_001365196', 'LOC645188', 'COL4A2-AS2', 'PPP5D1',
                    'FTCD-AS1', 'C9orf62', 'LOC107986453']
# free text no longer supported i.e. LRRK2:NM_198578.1:exon48:c.G7134C:p.K2378N
PROTEIN_FREE_TEXT_INITIALIZATION_RE = rf'([A-Z\d]*):(NM_[\d]*)[\.]?[\d]*:'

#  MUTATIONS CONSTANTS

MUTATION_REGEX = rf'p\.(?P<symbol>(?P<orig>[{A_A}]){{1}}(?P<location>[\d]+)(?P<change>[{A_A}]){{1}})'
REF_SEQ_PADDING = 5
NEW_MUTATION_DATA = {'chr': None, 'ref_na': None, 'alt_na': None, 'start': None, 'end': None, AFM_SCORE: NO_SCORE,
                     EVE_SCORE: NO_SCORE, ESM_SCORE: None, ESM_TYPE: NO_TYPE, ESM3_SCORE: None, ESM3_TYPE: NO_TYPE,
                     EVE_PREDICTION: NO_SCORE, EVE_TYPE: NO_TYPE, DS_RANK: None}
#  PATIENTS AND FAMILY

DEFAULT_PATIENT_COLUMNS = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Protein', 'Variant']

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
EVE_SCORE_COLUMN = "EVE_scores_ASM"
EVE_PREDICTION_COLUMN = "EVE_classes_75_pct_retained_ASM"
EVE_MUTATION_COLUMN = 'mt_aa'

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
CPT_NO_EVE_DATA_PATH = pjoin(CPT_PATH, 'CPT1_score_no_EVE_set')
CPT_EVE_IMPUTE_PATH_1 = pjoin(CPT_PATH, CPT_IMPUTE_DATA_1_NAME, CPT_IMPUTE_DATA_1_NAME)
CPT_EVE_IMPUTE_PATH_2 = pjoin(CPT_PATH, CPT_IMPUTE_DATA_2_NAME, CPT_IMPUTE_DATA_2_NAME)
CPT_INGENE_PATH = pjoin(CPT_PATH, CPT_EVE_DATA_NAME, CPT_EVE_DATA_NAME)
CPT_EXGENE_PATH = pjoin(CPT_PATH, CPT_IMPUTE_DATA_NAME, CPT_IMPUTE_DATA_1_NAME)
CPT_SCORE_COLUMN = "CPT1_score"
CPT_MUTATION_COLUMN = 'mutant'
CPT_TRUNCATE_TAIL = 7

#  ESM

ESM_AA_ORDER = 'LAGVSERTIDPKQNFYMHWC'
ESM_AA_EXTENDED = 'LAGVSERTIDPKQNFYMHWCXBUZO.-|'
AA_ESM_LOC = {aa: idx for idx, aa in enumerate(ESM_AA_ORDER)}
ESM_MAX_LENGTH = 1020

#  ESM-1b

#ESM1_AA_ORDER = 'KRHEDNQTSCGAVLIMPYFW'
#AA_ESM1_LOC = {aa: idx for idx, aa in enumerate(ESM1_AA_ORDER)}
ESM1B_MODEL = 'esm1b_t33_650M_UR50S'
ESM_VARIANTS_DATA = 'https://huggingface.co/spaces/ntranoslab/esm_variants/resolve/main/ALL_hum_isoforms_ESM1b_LLR.zip'
ESM_DATA = 'ESM_1b_variants'
ESM_DATA_PATH = pjoin(ESM_PATH, ESM_DATA)
ESM_INDEX_PATH = pjoin(ESM_PATH, 'index.json')
ESM_VARIANTS_PATH = pjoin(ESM_DATA_PATH, 'content', 'ALL_hum_isoforms_ESM1b_LLR')
ESM_PARQUET_DIR = Path(ESM_PATH) / "parquet_files"

ESM_FILE_SUFFIX = '_LLR.csv'
MASK_TOKEN = '<mask>'
REP_LAYERS = [33]

# ESM3

ESM3_MODEL = 'esm3_sm_open_v1'

#  DSRANK AND SUMMARY
EVE_COL, EVE_TYPE_COL, ESM_COL, ESM_TYPE_COL, ESM3_COL, ESM3_TYPE_COL, AFM_COL, DS_COL = \
    "eve", "eve_type", "esm1b", "esm1b_type", "esm3", "esm3_type", "afm", "ds_rank"
PROT_COL, MUT_COL = 'protein', 'variant'
COLUMNS_W_STATUS = [PROT_COL, MUT_COL, EVE_COL, EVE_TYPE_COL, ESM_COL, ESM_TYPE_COL, ESM3_COL, ESM3_TYPE_COL,
                    AFM_COL, DS_COL]
COLUMNS_NO_STATUS = [PROT_COL, MUT_COL, EVE_COL, ESM_COL, ESM3_COL, AFM_COL, DS_COL]
