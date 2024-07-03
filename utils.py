import os
from pathlib import Path
import glob
import warnings
import requests
from requests.adapters import HTTPAdapter, Retry
from definitions import *
import sys
import psutil
import pandas as pd
import gzip
import shutil
import hashlib
from tqdm import tqdm
import click


ALPHAFOLD_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v1.pdb"


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


def make_fasta(path, name, seq):
    full_path = os.path.join(path, f"{name}.fasta")
    with open(full_path, "w+") as file:
        file.write(f">{name}\n")
        file.write(f"{seq}\n")


def adaptive_chunksize(rowsize, ram_usage=0.5):
    """
    
    :param rowsize: float size of dataframe row in bits
    :param ram_usage: float [0,1] portion of available ram to use.
    Note setting ram_usage = 1.0 may result in memory errors.
    :return: int number of rows in chunk
    """
    available = psutil.virtual_memory().available * ram_usage
    return available // rowsize


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


def generate_proteins(proteins):
    """
    yields Protein from list of protein names =
    """
    for p_name in proteins:
        yield Protein.Protein(ref_name=p_name)


def generate_pairs(pairs):
    """
    yields Mutations from list of tuples ((protein_name, mutation_name),(protein_name, mutation_name))
    """
    for m1, m2 in pairs:
        m1 = Mutation.Mutation(m1[1], m1[0])
        m2 = Mutation.Mutation(m2[1], m2[0])
        yield Pair.Pair(m1, m2)


def generate_mutations(mutations):
    """
    yields Mutations from list of tuples [(mutation_name, protein_name)]
    """
    for m_name, p_name in mutations:
        yield Mutation.Mutation(m_name, p_name)


def create_mutation(s):
    s = s.split("_")
    return Mutation.Mutation(f"p.{s[1]}", s[0])


def update_patients_ranks():
    for p in [os.path.basename(i)[:-4] for i in glob.glob('DB/Patients/HR*')]:
        print(f"processing {p}")
        patient = Patient.Patient(p)
        f_id = p.split('_')[0]
        family = Family.Family(family_id=f_id)
        rank_1 = patient.top_all()
        rank_2 = family.recurring
        rank_3 = patient.in_n_params(2)
        rank_4 = (set(patient.top_interface()) & patient.in_n_params(1, interface=0)) | \
                 patient.in_n_params(3) | (family.recurring & patient.in_n_params(2)) | family.multiple_recurrence()
        ranks = []
        for mut in patient.mutations():
            if mut in rank_4:
                ranks.append(4)
                continue
            if mut in rank_3:
                ranks.append(3)
                continue
            if mut in rank_2:
                ranks.append(2)
                continue
            if mut in rank_1:
                ranks.append(1)
                continue
            else:
                ranks.append(0)
        patient._data['rank'] = ranks
        patient._data.to_csv(rf'DB/Patient/new/{p}')


def download_af_pdb(id):
    pass


def protein_exists(ref_name):
    return ref_name in set(os.listdir('DB/proteins'))


def all_but_one(patients):
    """
    input: list of patients
    output: set of mutation in n-1 patients
    """
    final_set = set()
    if len(patients) == 1:
        return set(patients[0].mutations())
    if len(patients) == 2:
        return set(patients[0].mutations()) | set(patients[1].mutations())
    for i in range(len(patients)):
        cur_patients = patients[0:i] + patients[i + 1:len(patients)]
        sets = [set(p.mutations()) for p in cur_patients]
        final_set |= set.intersection(*sets)
    return final_set


def recurring_family_mutations(
        families=['HR1', 'HR3', 'HR4', 'HR5', 'HR6', 'HR7', 'HR8', 'HR9', 'HR10', 'HR11', 'HR12']):
    mult_rec = set()
    for i in range(len(families)):
        f1 = Family.Family(families[i])
        for id2 in families[i + 1:]:
            f2 = Family.Family(id2)
            mult_rec |= f1.recurring & f2.recurring
    return mult_rec


def add_record(prot):
    if not prot.mutations:
        return
    for mut, details in prot.mutations.items():
        data['name'].append(prot.name)
        data['variant'].append(mut)
        data['eveScore'].append(details['eveScore'])
        data['evePrediction'].append(details['evePrediction'])
        data['firmScore'].append(details['firmScore'])
        if details['bertScore'] != -1:
            for i, score in enumerate(details['bertScore']):
                data[f'bert_{i + 1}'].append(score)
            data['bert_sum'].append(sum(list(details['bertScore'])))
            data['bert_max'].append(max(list(details['bertScore'])))
        else:
            for i in range(1, 6):
                data[f'bert_{i}'].append(-1)
            data['bert_sum'].append(-1)
            data['bert_max'].append(-1)
        data['mentioned'].append(str(Analyze.ProteinAnalyzer.search_litriture(prot.name)))


def find_mutations_with_pdbs(protein, log="log.txt", found='found.txt'):
    """
    iterates ove protein mutations to find those with pdbs
    :param protein: protein object
    :return: set (mutation names)
    """
    for mutation in protein.mutations:
        mut_obj = Mutation.Mutation(mutation, protein)
        if mut_obj.has_pdbs():
            print(f"{protein.name} - {mut_obj.name} have pdbs {mut_obj.pdbs}")
            with open(found, 'a') as file:
                file.write(f"{protein.name} - {mut_obj.name} - {mut_obj.pdbs}\n")
        else:
            with open(log, 'a') as file:
                file.write(f"{protein.name} - {mut_obj} - {mut_obj.pdbs}\n")


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

'''
def _firm_setup(ref_gen='GRCh37', n_threads=4):
    """
    sets up firm network (might take up to 15 min)
    :param ref_gen: reference genome to work on
    :param n_threads: number of threads to runt the process on
    :return: firm_classifier, geneffect_setup, thread_pool
    """
    thread_pool = multiprocessing.Pool(n_threads)
    geneffect_setup = geneffect.Setup(ref_gen)  # Must be a human reference genome
    firm.setup_uniprot_tracks(geneffect_setup)
    firm_classifier = firm.load_classifier(geneffect_setup)
    return firm_classifier, geneffect_setup, thread_pool
'''
