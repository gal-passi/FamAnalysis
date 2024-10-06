import glob
from Bio import Entrez, PDB
from urllib.error import HTTPError as HTTPError
from definitions import *
from utils import print_if, warn_if, create_session, safe_post_request, safe_get_request
from utils import progress_bar, process_fastas, read_in_chunks, process_pdb_ids_query
import wget
from os.path import basename


class Uniport:
    """
    This class is responsible to connect to online DBs and retrieve information
    """

    def __init__(self, verbose_level=1):
        Entrez.email = CONTACT
        self._v = verbose_level
        self.pdpl = PDB.PDBList()

    def fetch_uniport_sequences(self, uid, ):
        """
        Retrieve all known isoforms from uniport
        :param uid: Uniport id only primary name
        :return: {uid_iso_index: sequence}
        """
        print_if(self._v, VERBOSE['thread_progress'], "Retrieving isoforms from Uniprt...")
        s = create_session(DEFAULT_HEADER, RETRIES, WAIT_TIME, RETRY_STATUS_LIST)
        url = Q_UNIP_ALL_ISOFORMS.format(uid, uid)
        response = safe_post_request(s, url, TIMEOUT, self._v, CON_ERR_FUS.format(uid, url))
        if not response.ok:
            print_if(self._v, VERBOSE['thread_warnings'], CON_ERR_GENERAL.format('fetch_uniprot_sequences', uid))
            return {}
        if response.text == '':
            print_if(self._v, VERBOSE['thread_warnings'], f"no sequences found for {uid}")
            return {}
        return process_fastas(response.text)

    def expend_isoforms(self, prot, limit=20, unique_key=''):
        """
        will add all known isoforms of prot.name including non-reviewed
        will be returned in form {uid_iso:...}
        this will noe override the default protein isoforms
        :param prot: protein obj
        :param limit: max number of uids to quary
        :param unique_key: str added to search problematic proteins
        :return: {uid_iso_index: seq}
        """
        isoforms = {}
        if not unique_key:
            query = UNIPORT_QUERY_URL + Q_ISOFORMS_PROT.format(prot.name)
        else:
            query = UNIPORT_QUERY_URL + Q_ISOFORMS_KEY.format(unique_key)

        s = create_session(DEFAULT_HEADER, RETRIES, WAIT_TIME, RETRY_STATUS_LIST)
        r = safe_get_request(s, query, TIMEOUT, self._v, CON_ERR_EI.format(prot.name))
        if not r:
            return {}
        if r.text == '':
            return {}
        uids = r.text.split("\n")
        print_if(self._v, VERBOSE['thread_progress'], f"found {len(uids) - 2} uids")
        for uid in r.text.split("\n")[1:limit + 1]:
            if uid == '':
                return isoforms
            res = self.fetch_uniport_sequences(uid, expend=True)
            isoforms = {**isoforms, **res}
        return isoforms

    def fetch_pdbs(self, uids=None, prot=None, reviewed=True):
        """
        retreives all known BPDs of a the given uniport id or protein object
        :param uid: set of uniprot ids
        :param prot: optional Protein obj find pdbs for all known uids
        :param reviewed: bool, default True, if False non-reviewed pdbs will be added.
                NOTE this may exyended running time significantly.
        :return: dict {pdb_id: seqence}
        """
        if prot:
            reviewed_ids = set(prot.all_uids()['reviewed'])
            nreviewed_ids = set(prot.all_uids()['non_reviewed']).difference(reviewed_ids)
            uids = reviewed_ids if reviewed else nreviewed_ids
        elif uids:
            if isinstance(uids, str):
                uids = {uids}
        else:
            raise ValueError(f"failed on {uids}")
        print_if(self._v, VERBOSE['thread_progress'], f"Fetching Pdb ids...")
        return self._fetch_pdbs_in_chunks(list(uids))


    def _fetch_pdbs_in_chunks(self, uids):
        """
        helper for fetch pdbs
        :param uid: set of uniprot ids to fetch
        :return: dict {pdb_id: sequence}
        """
        pdbs = {}
        s = create_session(DEFAULT_HEADER, RETRIES, WAIT_TIME, RETRY_STATUS_LIST)
        for chunk in read_in_chunks(uids, CHUNK_SIZE):
            query = Q_UNIP_ALL_PDBS
            format = ''.join([f"((accession:{id})+OR+)" for id in chunk])
            if not format:
                return pdbs
            format = format[:-4]  # remove last OR
            query = query.format(format)
            r = safe_get_request(s, query, TIMEOUT, self._v, CON_ERR_FP_1.format(query))
            if not r.text:
                continue
            pdb_ids = process_pdb_ids_query(r.text)
            pdbs |= self._query_pdb_seqs(s, pdb_ids)

    @staticmethod
    def _query_pdb_seqs(session, ids):
        """
        helper for fetch pdbs
        :param ids: list of pdb ids to fetch
        :return: dict {pdb_id: sequence}
        """
        delta = len(ids) / 2  # for proteins with many pdbs default timeout might not be enough
        pdbs = {}
        for id in ids:
            url = EBI_PDB_URL + f'{id}'
            value = safe_get_request(session, url, TIMEOUT + delta, CON_ERR_FP_2.format(id, id))
            if not value:
                continue
            pdbs[id] = value.json()[id.lower()][0]['pdb_sequence']
        return pdbs


    def uid_from_name(self, name):
        """
        return uniprot-id given a protein ref_name
        :param name: protein ref_name
        :param all: optional return a list of all uids found
        :param reviewed: default is True will return only reviewed entities. set to False if no results are found
        :return: {'reviewed': [...], 'non_reviewed': [...]}
        """
        ret = {'reviewed': [], 'non_reviewed': [], 'main_entery': [], 'all_enteries': [], 'aliases': []}
        query = UNIPORT_QUERY_URL + Q_UID_PROT_ALL.format(name)
        s = create_session(DEFAULT_HEADER, RETRIES, WAIT_TIME, RETRY_STATUS_LIST)
        # r = req.get(query, timeout=TIMEOUT)
        r = safe_get_request(s, query, TIMEOUT, self._v, CON_ERR_UFN.format(name))
        if not r:
            return ret
        if r.text == '':
            return ret
        ret = self._process_uid_query(r.text)
        return ret

    @staticmethod
    def _process_uid_query(data):
        ret = {'reviewed': [], 'non_reviewed': [], 'main_entery': [], 'all_enteries': [], 'aliases': []}
        rows = data.split('\n')[1:-1]  # first is header last is blank
        if not rows:
            return ret
        main_entry = rows[0].split('\t')[UIDS_COL_IDX]  # first entry is considered main
        ret['main_entery'].append(main_entry)
        for row in rows:
            values = row.split('\t')
            entry, reviewed, gene = values[UIDS_COL_IDX], values[REVIEWED_COL_IDX], values[GENE_NAME_COL_IDX].split(" ")
            ret['all_enteries'].append(entry)
            ret['aliases'] += gene
            if reviewed == UNIP_REVIEWED:
                ret['reviewed'].append(entry)
            else:
                ret['non_reviewed'].append(entry)
        return ret

    def alphafold_confidence(self, prot, mut):
        """
        returns confidence of alphafold model - if unable to find or model or sequences don't match return -1
        """
        s = create_session(DEFAULT_HEADER, RETRIES, WAIT_TIME, RETRY_STATUS_LIST)
        # resp = req.get(ALPHAFOLD_PDB_URL.format(prot.Uid), timeout=TIMEOUT)
        resp = safe_get_request(s, ALPHAFOLD_PDB_URL.format(prot.Uid), TIMEOUT, self._v,
                                CON_ERR_GENERAL.format('alphafold_confidence', prot.Uid))
        if not resp:
            return -1
        if not resp.ok:
            warn_if(self._v, VERBOSE['thread_warnings'], f"Failed to find alphafold model for {prot.Uid}")
            return -1
        relavent_confidence = [float(i[61:66]) for i in resp.content.split(b"\n") if
                               i.startswith(b"ATOM  ") and int(i[22:26]) == mut.loc and
                               AA_SYN[mut.origAA] == i[16:20].decode('utf-8').strip()]

        if len(relavent_confidence) < 1:
            # try to find in sequence assumes reference length of 10
            if not mut.ref_seqs:
                warn_if(self._v, VERBOSE['thread_warnings'],
                        f"Failed to find residue {mut.origAA} in {mut.loc} -- no ref seqs")
                return -1
            sequence = self._obtain_seq(resp)
            ref = next(iter(mut.ref_seqs.values()))
            index = sequence.find(ref)
            if index == -1:
                warn_if(self._v, VERBOSE['thread_warnings'],
                        f"Faild to find residue {mut.origAA} in {mut.loc} -- can't find reference in sequence")
                return -1
            # if reference sequence is not 10 need to change "5"
            relavent_confidence = [float(i[61:66]) for i in resp.content.split(b"\n") if
                                   i.startswith(b"ATOM  ") and int(i[22:26]) == index + 6]

        return relavent_confidence[0]

    def alpha_seq(self, prot):
        s = create_session(DEFAULT_HEADER, RETRIES, WAIT_TIME, RETRY_STATUS_LIST)
        resp = safe_get_request(s, ALPHAFOLD_PDB_URL.format(prot.Uid), TIMEOUT, self._v,
                                CON_ERR_GENERAL.format('alpha_seq', prot.Uid))
        if not resp:
            return {}
        if not resp.ok:
            return {}
        return {'alpha': self._obtain_seq(resp)}

    def _obtain_seq(self, bytefile):
        """
        obtaipns protein sequence from bytefile response
        """
        seqs = [row[17:].split(b' ')[2:15] for row in bytefile.content.split(b"\n") if row.startswith(b"SEQRES ")]
        seqs = [AA_SYN_REV[item.decode('utf-8')] for sublist in seqs for item in sublist if item != b'']
        return "".join(seqs)

    def download_pdb(self, id, path):
        """
        downloads specific pdb file and save to path as id.pdb
        :param id: str id of pdb file/cs/usr/gal_passi/dina-lab/sourc
        :param path: str directory to save the file to
        :return:
        """
        self.pdpl.retrieve_pdb_file(id, pdir=path, file_format="pdb")

    def download_eve_data(self, prot_name):
        """
        :param prot_names: str
        :return: bool True is successful
        """
        s = create_session(DEFAULT_HEADER, RETRIES, WAIT_TIME, RETRY_STATUS_LIST)
        existing_files = {basename(p)[:-4] for p in glob.glob(pjoin(EVE_VARIANTS_PATH, '*.csv'))}
        name = prot_name + '_HUMAN'
        if name in existing_files:
            return True
        url = EVE_SINGLE_PROTEIN.format(name)
        if s.get(url, timeout=TIMEOUT).ok:
            print_if(self._v, VERBOSE['thread_progress'], EVE_PROT_DOWNLOAD_MSG.format(prot_name))
            for _ in RETRIES:
                try:
                    wget.download(url, bar=progress_bar, out=EVE_DATA_PATH)
                    return True
                except:
                    pass
            return False
        else:
            return False


# will be deprecated

def fatch_all_NCBIs(self, ncbi_id):
    idx = 0
    united = {}
    while True:
        record = self.fetch_NCBI_seq(ncbi_id + f".{idx}")
        if (len(record) == 0) and (idx > 0):
            return united
        united = {**united, **record}
        idx += 1


def fetch_NCBI_seq(self, ncbi_id):
    """
    fetches amino acid sequence from ncbi
    :param ncbi_id: str (format NM_001282166.1)
    :return: {ncbi_id : sequence}
    """
    print_if(self._v, VERBOSE['thread_progress'], f"Retrieving isoform {ncbi_id} from NCBI...")
    try:
        handle = Entrez.efetch(db="nucleotide", id=ncbi_id, retmode="xml")
    except HTTPError as e:
        warn_if(self._v, VERBOSE['thread_warnings'],
                f'Unable to fetch NCBI record with id {ncbi_id} HTTPError...skipping')
        warn_if(self._v, VERBOSE['raw_warnings'], f'\n{e}')
        return {}
    records = Entrez.read(handle)
    try:
        for dict_entry in records[0]["GBSeq_feature-table"]:
            if dict_entry['GBFeature_key'] == 'CDS':
                for sub_dict in dict_entry["GBFeature_quals"]:
                    if sub_dict['GBQualifier_name'] == 'translation':
                        print_if(self._v, VERBOSE['thread_progress'], f"done")
                        return {ncbi_id: sub_dict["GBQualifier_value"]}
    except Exception as e:
        warn_if(self._v, VERBOSE['thread_warnings'], f'Unable to fetch NCBI record with id {ncbi_id}...skipping')
        warn_if(self._v, VERBOSE['raw_warnings'], f'\n{e}')
        return {}
