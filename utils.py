import os
from pathlib import Path
import warnings
import requests
import torch.nn
from requests.adapters import HTTPAdapter, Retry
from definitions import *
import sys
import psutil
import pandas as pd
import gzip
import hashlib
from tqdm import tqdm
import click
import re
import tempfile
from Bio import SeqIO
#from esm.models.esm3 import ESM3
#from esm.pretrained import ESM3_sm_open_v0
#from esm.sdk.api import ESM3InferenceClient
#from esm.tokenization.sequence_tokenizer import EsmSequenceTokenizer
#from huggingface_hub import login



def print_if(verbose: object, thr: object, text: object) -> object:
    """
    print text if verbose > thr
    :param verbose: int
    :param thr: int
    :param text: str
    :return:
    """
    if verbose >= thr:
        print(text)


def warn_if(verbose, thr, text):
    """
    print text if verbose > thr
    :param verbose: int
    :param thr: int
    :param text: str
    :return:
    """
    if verbose >= thr:
        warnings.warn(text)


def create_session(header, retries=5, wait_time=0.5, status_forcelist=None):
    """
    Creates a session using pagination
    :param header: str url header session eill apply to
    :param retries: int number of retries on failure
    :param wait_time: float time (sec) between attempts
    :param status_forcelist: list HTTP status codes that we should force a retry on
    :return: requests session
    """
    s = requests.Session()
    retries = Retry(total=retries,
                    backoff_factor=wait_time,
                    status_forcelist=status_forcelist)

    s.mount(header, HTTPAdapter(max_retries=retries))
    return s


def progress_bar(current, total, width=80):
    progress_message = "Downloading: %d%% [%d / %d] bytes" % (current / total * 100, current, total)
    sys.stdout.write("\r" + progress_message)
    sys.stdout.flush()


def safe_get_request(session, url, timeout, verbose_level, warning_msg='connection failed', return_on_failure=None,
                     warning_thr=VERBOSE['thread_warnings'], raw_err_thr=VERBOSE['raw_warnings']):
    """
    creates a user friendly request raises warning on ConnectionError but will not crush
    verbose_level = 3 will return raw Error massage in warning
    :param raw_err_thr: int threshold to print raw error messages
    :param warning_thr: int threshold to print warning messages
    :param session: requests session obj
    :param url: str url to query
    :param timeout: float max time to wait for response
    :param verbose_level: int
    :param warning_msg: str msg to display on failure
    :param return_on_failure: value to return upon exception
    :return: response
    """
    try:
        r = session.get(url, timeout=timeout)
    except requests.exceptions.ConnectionError as e:
        warn_if(verbose_level, warning_thr, warning_msg)
        warn_if(verbose_level, raw_err_thr, f"{e}")
        return return_on_failure
    return r


def safe_post_request(session, url, timeout, verbose_level, warning_msg='connection failed', return_on_failure=None,
                      warning_thr=VERBOSE['thread_warnings'], raw_err_thr=VERBOSE['raw_warnings']):
    """
    creates a user friendly request raises warning on ConnectionError but will not crush
    verbose_level = 3 will return raw Error massage in warning
    :param session: requests session obj
    :param url: str url to query
    :param timeout: float max time to wait for response
    :param verbose_level: int
    :param warning_msg: str msg to display on failure
    :param return_on_failure: value to return upon exception
    :param raw_err_thr: int threshold to print raw error messages
    :param warning_thr: int threshold to print warning messages
    :return: response
    """
    try:
        r = session.post(url, timeout=timeout)
    except requests.exceptions.ConnectionError as e:
        warn_if(verbose_level, warning_thr, warning_msg)
        warn_if(verbose_level, raw_err_thr, f"{e}")
        return return_on_failure
    return r


def read_and_check_response(resp, query):
    if not resp:
        raise ConnectionError(f'querry: {query} --> Failed')
    if not resp.ok:
        raise ConnectionError(f'querry: {query} --> Failed with code {resp.status_code}')
    return resp.text


def remove_whitespaces(text):
    """
    removes all whitespaces from string
    :param str:
    :return:
    """
    return re.sub('\s', '', text)


def process_fastas(text):
    """
    process multiple fasta sequences
    :return: {id: sequence}
    """
    temp = tempfile.TemporaryFile(mode='w+t')
    temp.writelines(text)
    temp.seek(0)
    ret = {seq_record.id: str(seq_record.seq) for seq_record in SeqIO.parse(temp, "fasta")}
    temp.close()
    return ret


def process_pdb_ids_query(text):
    ids = []
    for row in text.splitlines()[1:]:
        ids += row.split('\t')[-1].split(';')[:-1]
    return ids


def make_fasta(path, name, seq):
    full_path = os.path.join(path, f"{name}.fasta")
    with open(full_path, "w+") as file:
        file.write(f">{name}\n")
        file.write(f"{seq}\n")


def na_to_aa_translator(na_seq):
    assert len(na_seq) % 3 == 0, 'na seq must be a multiple of 3'
    aa_seq = ''
    for idx in range(0, len(na_seq), 3):
        codon = na_seq[idx: idx + CODON_LENGTH].lower()
        if codon in STOP_CODONS:
            print(f'stop codon at {idx} out of {len(na_seq)} --> {codon}')
            break
        aa_seq += CODON_TRANSLATOR[codon]
    print(f'na_length: {len(na_seq) // 3} \t aa_length: {len(aa_seq)}')
    assert (len(na_seq) // 3) - 1 == len(aa_seq), 'na seq and aa seq length do not match - check for early stop codon'
    return aa_seq


def adaptive_chunksize(rowsize, ram_usage=0.5):
    """
    
    :param rowsize: float size of dataframe row in bits
    :param ram_usage: float [0,1] portion of available ram to use.
    Note setting ram_usage = 1.0 may result in memory errors.
    :return: int number of rows in chunk
    """
    available = psutil.virtual_memory().available * ram_usage
    return available // rowsize


def read_in_chunks(array, chunk_size):
    for i in range(0, len(array), chunk_size):
        yield array[i:i + chunk_size]


def afm_iterator(chunksize, usecols=None):
    """
    :param chunksize: int
    :param usecols: Sequence of Hashable or Callable, optional. Subset of columns to select,
    :yields: DataFrame of size chuksize
    """
    for chunk in pd.read_csv(AFM_DATA_PATH, sep='\t', chunksize=chunksize,
                             header=AFM_HEADER, usecols=usecols):
        yield chunk


def afm_range_read(idx_from, idx_to, usecols=None):
    nrows = idx_to - idx_from
    return pd.read_csv(AFM_DATA_PATH, sep='\t', names=AFM_COL_NAMES, header=AFM_HEADER, skiprows=idx_from,
                       nrows=nrows, usecols=usecols)


def ugzip(path, outfile, chunksize):
    """
    unzips .gz file in chunks
    :param path: str path to .gz file
    :param outfile: str outfile path
    :return:
    """
    chunksize = int(chunksize)
    with gzip.open(path, 'rb') as f_in:
        with open(outfile, 'wb') as f_out:
            chunk = f_in.read(chunksize)
            while chunk:
                f_out.write(chunk)
                chunk = f_in.read(chunksize)


def name_for_esm(name):
    if '_HUMAN' in name:
        return name[:-6]
    return name


def name_for_cpt(name):
    if name.endswith('_HUMAN'):
        return name
    return name + '_HUMAN'


def sequence_from_esm_df(esm_data):
    return ''.join([desc[0] for desc in esm_data.columns.to_list()[1:]])


def sequence_from_cpt_df(cpt_data):
    return ''.join([aa_change[0] for aa_change in cpt_data['mutant'][0::20]])


def protein_exists(ref_name):
    return ref_name in set(os.listdir('DB/proteins'))


def seacrh_by_ref_seq(main_seq, ref_seq, padding=REF_SEQ_PADDING):
    """
    :return: index of start row of wt AA
    """
    seq_start = main_seq.find(ref_seq)  # TODO handle edges - will be implemented in Mutation
    if seq_start == -1:
        return -1
    aa_location = seq_start + padding
    start_row_idx = aa_location * 20
    return start_row_idx


def summary_df(include_status=False):
    """
    :param include_status: bool include esm and eve score type
    :return: DataFrame template for summary
    """
    if include_status:
        return pd.DataFrame(columns=COLUMNS_W_STATUS)
    else:
        return pd.DataFrame(columns=COLUMNS_NO_STATUS)

def esm_setup(model_name=ESM1B_MODEL):
    """
    :param model_name: str model name
    :return: model, alphabet api
    """
    print(model_name)
    #model, alphabet = torch.hub.load("facebookresearch/esm:main", 'esm1b_t33_650M_UR50S', source="github", force_reload=True)
    model, alphabet = torch.hub.load("facebookresearch/esm:main", 'esm1b_t33_650M_UR50S')
    model = model.to(DEVICE)
    print(f"model loaded on {DEVICE}")
    return model, alphabet

def esm3_setup(model_name=ESM3_MODEL, token=HUGGINGFACE_TOKEN):
    """
    :param model_name: str model name
    :param token: huggingface login token
    :return: ESM3_model, ESM_tokenizer
    """
    login(token)
    model: ESM3InferenceClient = ESM3.from_pretrained(model_name).to(DEVICE)
    if DEVICE == 'cuda':
        model = model.float()
    print(f"model loaded on {DEVICE}")
    tokenizer = EsmSequenceTokenizer()
    return model, tokenizer

def esm_seq_logits(model, tokens, log=True, softmax=True, return_device='cpu', esm_version=3):
    """
    :param model: ESM3 initiated model
    :param tokens: tokenized sequence
    :param log: bool return log of logits
    :param softmax: bool return softmax of logits
    :param return_device: 'cpu' | 'cuda should logits be returned on cpu or cuda
    :return: numpy array - log_s
    """

    if return_device == 'cuda':
        assert return_device == DEVICE, 'return type cuda but no GPU detected'
    model.eval()
    with torch.no_grad():
        if esm_version == 3:
            assert tokens.device == model.device, 'tokens and model must be on the same device'
            output = model.forward(sequence_tokens=tokens)
            sequence_logits = output.sequence_logits
        else:  #  esm 1,2
            assert tokens.device == next(model.parameters()).device, 'tokens and model must be on the same device'
            sequence_logits = model(tokens, repr_layers=REP_LAYERS, return_contacts=False)['logits']
        aa_logits = sequence_logits[0, 1:-1, 4:24]
        if log and softmax:
            token_probs = torch.log_softmax(aa_logits, dim=-1)
        elif log:
            token_probs = torch.log(aa_logits)
        elif softmax:
            token_probs = torch.softmax(aa_logits, dim=-1)
        token_probs = token_probs.cpu().numpy()
        return token_probs


def esm_process_long_sequences(seq, loc):
    """
    trims seq to len < 1024 using mut Object. to fit for bert model
    :return: offset from orig location, trimmed sequence
    """
    thr = ESM_MAX_LENGTH // 2
    left_bound = 0 if loc-thr < 0 else loc - thr
    right_bound = len(seq) if loc + thr >len(seq) else loc+thr
    left_excess = 0 if loc - thr > 0 else abs(loc-thr)
    right_excess = 0 if loc + thr <= len(seq) else loc+thr-len(seq)
    if (left_excess == 0) and (right_excess == 0):
        return left_bound, seq[left_bound:right_bound]
    if left_excess > 0:
        return 0, seq[left_bound:right_bound+left_excess]
    if right_excess > 0:
        return left_bound-right_excess, seq[left_bound-right_excess:right_bound]

def create_patient_data(data, patient, columns, empty=None):
    """
    dataframe
    :param columns: list of columns to include in patient DataFrame
    :param data: DataFrame containing all data
    :param patient: patient column name
    :param empty: value for empty value if None assume Na
    :return: df patient data
    """
    return data[data[patient] is not empty][columns]


class SafeDownloader:
    """
    code by tobiasraabe - Tobias Raabe
    cloned from https://gist.github.com/58adee67de619ce621464c1a6511d7d9.git
    """

    def __init__(self, urls, file_names, url_hashes=None, outfile='.', url_base='', block_size=1024, verbose_level=1):
        """
        :param urls: list of urls to download,
        :param url_hashes: hash or urls in format sha256 lowercase
        :param file_names: list of strings names save names of files
        :param outfile: string directory to download to
        :param url_base: optional - shared url base for readability
        :param block_size: int size in bits of download blocks
        """
        self.base = url_base
        self.urls = urls
        self.url_hashes = url_hashes
        self.file_names = file_names
        self.outfile = Path(outfile)
        self.block_size = block_size
        self._v = verbose_level
        self._context_setting = dict(help_option_names=['-h', '--help'])

    def downloader(self, position: int, resume_byte_pos: int = None):
        """Download url in ``URLS[position]`` to disk with possible resumption.

        Parameters
        ----------
        position: int
            Position of url.
        resume_byte_pos: int
            Position of byte from where to resume the download

        """
        # Get size of file
        url = self.urls[position]
        r = requests.head(url)
        file_size = int(r.headers.get('content-length', 0))

        # Append information to resume download at specific byte position
        # to header
        resume_header = ({'Range': f'bytes={resume_byte_pos}-'}
                         if resume_byte_pos else None)

        # Establish connection
        r = requests.get(url, stream=True, headers=resume_header)

        # Set configuration
        initial_pos = resume_byte_pos if resume_byte_pos else 0
        mode = 'ab' if resume_byte_pos else 'wb'
        file = self.outfile / self.file_names[position]

        with open(file, mode) as f:
            with tqdm(total=file_size, unit='B',
                      unit_scale=True, unit_divisor=self.block_size,
                      desc=file.name, initial=initial_pos,
                      ascii=True, miniters=1) as pbar:
                for chunk in r.iter_content(32 * self.block_size):
                    f.write(chunk)
                    pbar.update(len(chunk))

    def download_file(self, position: int) -> None:
        """Execute the correct download operation.

        Depending on the size of the file online and offline, resume the
        download if the file offline is smaller than online.

        Parameters
        ----------
        position: int
            Position of url.

        """
        # Establish connection to header of file
        url = self.urls[position]
        r = requests.head(url)

        # Get filesize of online and offline file
        file_size_online = int(r.headers.get('content-length', 0))
        file = self.outfile / self.file_names[position]

        if file.exists():
            file_size_offline = file.stat().st_size

            if file_size_online != file_size_offline:
                print_if(self._v, VERBOSE['program_warning'], DOWNLOAD_INCOMPLETE_WRN.format(file))
                self.downloader(position, file_size_offline)
            else:
                print_if(self._v, VERBOSE['program_progress'], DOWNLOAD_COMPLETE_MSG.format(file))
                pass
        else:
            print_if(self._v, VERBOSE['program_progress'], DOWNLOAD_START_MSG.format(file))
            self.downloader(position)

    def validate_file(self, position: int) -> None:
        """Validate a given file with its hash.

        The downloaded file is hashed and compared to a pre-registered
        has value to validate the download procedure.

        Parameters
        ----------
        position: int
            Position of url and hash.

        """
        file = self.outfile / self.file_names[position]
        try:
            hash = self.url_hashes[position]
        except IndexError:
            print_if(self._v, VERBOSE['program_warning'], DOWNLOAD_NO_HASH_ERR.format(file.name))
            return 0

        sha = hashlib.sha256()
        with open(file, 'rb') as f:
            while True:
                chunk = f.read(1000 * 1000)  # 1MB so that memory is not exhausted
                if not chunk:
                    break
                sha.update(chunk)
        try:
            assert sha.hexdigest() == hash
        except AssertionError:
            file = self.file_names[position]
            print_if(self._v, VERBOSE['program_warning'], DOWNLOAD_CORRUPTION_ERR.format(file))
        else:
            print_if(self._v, VERBOSE['program_progress'], DOWNLOAD_VALIDATION_MSG.format(file))

    @click.group(context_settings=CONTEXT_SETTINGS, chain=True)
    def cli(self):
        """Program for downloading and validating files.

        It is possible to run both operations consecutively with

        .. code-block:: shell

            $ python python-downloader.py download validate

        To download a file, add the link to ``URLS`` and its hash to ``HASHES`` if
        you want to validate downloaded files.

        """
        pass

    def download(self):
        """Download files specified in ``URLS``."""
        for position in range(len(self.urls)):
            self.download_file(position)

    def validate(self):
        """Validate downloads with hashes in ``HASHES``."""
        for position in range(len(self.urls)):
            self.validate_file(position)



