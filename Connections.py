import requests as req
from Bio import Entrez, PDB, SeqIO
from urllib.error import HTTPError as HTTPError
import re
from definitions import *
from utils import print_if, warn_if



class Uniport:
    """
    This class is responsible to connect to online DBs and retrieve information
    """
    HEADERS = {'User-Agent': 'Python {}'.format(CONTACT)}


    def __init__(self, verbose_level=1):
        Entrez.email = CONTACT
        self._v = verbose_level
        self.pdpl = PDB.PDBList()

    def fetch_uniport_sequences(self, Uid, expend=False):
        """
        Retrieve all known isoforms from uniport
        :param Uid: Uniport id only primary name
        :param expend: bool optional whether to add uid to sequence name
        :return: {uid_iso_index: sequence}
        """
        print_if(self._v, 3, "Retrieving isoforms from Uniprt...")
        if not Uid: return ''
        sequences = {}
        index = 1
        while True:
            currentUrl = f"{UNIPORT_URL}{Uid}-{index}.fasta"
            response = req.post(currentUrl, timeout=TIMEOUT)
            seq = response.text
            if not response.ok:
                print_if(self._v, 3, "done")
                return sequences
            if not expend:
                sequences[f"iso_{index}"] = self.remove_whitespaces(seq[seq.find("\n")+1:])  # remove header
            else:
                sequences[f"{Uid}_iso_{index}"] = self.remove_whitespaces(seq[seq.find("\n")+1:])
            index += 1

    def expend_isoforms(self, prot, limit=20, unique_key = ''):
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
            query = UNIPORT_QUERY_URL + \
                    f"fields=id&format=tsv&query=gene_exact:{prot.name}+AND+organism_id:9606"
        else:
            query = UNIPORT_QUERY_URL + \
                    f"fields=id&format=tsv&query={unique_key}+AND+organism_id:9606"

        r = req.get(query, timeout=TIMEOUT)
        if r.text == '':
            return {}
        uids = r.text.split("\n")
        print_if(self._v, 3, f"found {len(uids) - 2} uids")
        for uid in r.text.split("\n")[1:limit+1]:
            if uid == '':
                return isoforms
            res = self.fetch_uniport_sequences(uid, expend=True)
            isoforms = {**isoforms, **res}
        return isoforms

    def fetch_pdbs(self, Uid="", prot=None, reviewed=True):
        """
        retreives all known BPDs of a the given uniport id or protein object
        :param Uid: uniport id
        :param prot: optional Protein obj find pdbs for all known uids
        :param reviewed: bool, default True, if False non-reviewed pdbs will be added.
                NOTE this may exyended running time significantly.
        :return: dict {pdb_id: seqence}
        """
        if prot:
            pdbs = {}
            reviewed = set(prot.all_uids()['reviewed'])
            nreviewed = set(prot.all_uids()['non_reviewed']).difference(reviewed)
            for uid in reviewed:
                pdbs = {**pdbs, **self.fetch_pdbs(uid)}
            if not reviewed:
                for uid in nreviewed:
                    pdbs = {**pdbs, **self.fetch_pdbs(uid)}
            return pdbs

        if not Uid:
            return {}

        #params = {'format': 'tab', 'query': 'ID:{}'.format(Uid), 'columns': 'id,database(PDB)'}
        query = UNIPORT_QUERY_URL + \
                f"fields=id,xref_pdb&format=tsv&query={Uid}"
        print_if(self._v, 3, f"Fetching Pdb ids for {Uid}...")
        r = req.get(query, timeout=TIMEOUT)
        print_if(self._v, 3, f"done")
        pdbs = str(r.text).splitlines()[-1].split('\t')[-1].split(';')[:-1]
        print_if(self._v, 3, f"Fetching pdbs sequences...")
        try:
            delta = len(pdbs) / 2  # for proteins with many pdbs default timeout might not be enough
            ret = {id: req.get(EBI_PDB_URL + f'{id}', timeout=TIMEOUT + delta).json()[id.lower()][0]['sequence']
                   for id in pdbs}
            print_if(self._v, 3, f"done")
        except Exception as e:
            warn_if(self._v, 2, f"exception on fetch_uniport with {Uid} -> {pdbs}\nreturning empty...")
            warn_if(self._v, 3, f"exception on fetch_uniport with {Uid} -> {pdbs}\nreturning empty...\n{e}")
            ret = {}

        return ret

    def uid_from_name(self, name, all=False, reviewed=True):
        """
        return uniprot-id given a protein ref_name
        :param name: protein ref_name
        :param all: optional return a list of all uids found
        :param reviewed: default is True will return only reviewed entities. set to False if no results are found
        :return:
        """
        UIDS = 0
        if reviewed:
            query = UNIPORT_QUERY_URL + \
                    f"fields=&id&format=tsv&query={name}+AND+organism_id:9606+AND+reviewed:true"
            #params = {'format': 'tab', 'query': f"gene_exact:{name} AND organism:homo_sapiens AND reviewed:yes", 'columns': 'id'}
        else:
            query = UNIPORT_QUERY_URL + \
                    f"fields=&id&format=tsv&query={name}+AND+organism_id:9606+AND+reviewed:false"
            #params = {'format': 'tab', 'query': f"gene_exact:{name} AND organism:homo_sapiens", 'columns': 'id'}
        r = req.get(query, timeout=TIMEOUT)
        if r.text == '':
                return []
        rows = r.text.split('\n')
        uids = [row.split('\t')[UIDS] for row in rows[1:-1]]
        if len(uids) == 0:
            return ''
        return uids[0] if not all else uids

    def entery_name(self, protein=None, ref_name='', all_results=False):
        """

        :param protein: Protein obj
        :param all_results: returns list of all entery name found if False returns the first
        :return:
        """
        ENTERY_NAME = 1
        if not protein and not ref_name:
            raise ValueError("Must include eithe protein object or protein ref_name")

        name = ref_name if ref_name else protein.name
        url = f"https://rest.uniprot.org/uniprotkb/search?fields=&gene&format=tsv&query={name}+AND+organism_id:9606"
        r = req.get(url, timeout=TIMEOUT)
        if r.text == '':
            return ''
        rows = r.text.split('\n')
        entery_names = [row.split('\t')[ENTERY_NAME] for row in rows[1:-1]]
        if len(entery_names) == 0:
            return ''
        return entery_names[0] if not all_results else entery_names

    def synonms(self, protein=None, by_name=None):
        """
        :param name: protein object
        :return: list of all synonyms for protein name
        """
        if not protein and not by_name:
            return []
        search_term = protein.Uid if protein else by_name
        url = f"https://rest.uniprot.org/uniprotkb/search?fields=&gene&format=tsv&query={search_term}"
        r = req.get(url, timeout=TIMEOUT)
        if r.text == '':
            return []
        res = []
        for line in r.text.split('\n')[1:-1]:
            res += line.split('\t')[4].split(" ")
        return res


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
        print_if(self._v, 3, f"Retrieving isoform {ncbi_id} from NCBI...")
        try:
            handle = Entrez.efetch(db="nucleotide", id=ncbi_id, retmode="xml")
        except HTTPError as e:
            warn_if(self._v, 2, f'Unable to fetch NCBI record with id {ncbi_id} HTTPError...skipping')
            warn_if(self._v, 3, f'Unable to fetch NCBI record with id {ncbi_id} HTTPError...skipping\n{e}')
            return {}
        records = Entrez.read(handle)
        try:
            for dict_entry in records[0]["GBSeq_feature-table"]:
                if dict_entry['GBFeature_key'] == 'CDS':
                    for sub_dict in dict_entry["GBFeature_quals"]:
                        if sub_dict['GBQualifier_name'] == 'translation':
                            print_if(self._v, 3, f"done")
                            return{ncbi_id: sub_dict["GBQualifier_value"]}
        except Exception as e:
            warn_if(self._v, 2, f'Unable to fetch NCBI record with id {ncbi_id}...skipping')
            warn_if(self._v, 3, f'Unable to fetch NCBI record with id {ncbi_id}...skipping\n{e}')
            return {}

    def alphafold_confidence(self, prot, mut):
        """
        returns confidence of alphafold model - if unable to find or model or sequences don't match return -1
        """
        resp = req.get(ALPHAFOLD_PDB_URL.format(prot.Uid), timeout=TIMEOUT)
        if not resp.ok:
            warn_if(self._v, 2, f"Failed to find alphafold model for {prot.Uid}")
            return -1
        relavent_confidence = [float(i[61:66]) for i in resp.content.split(b"\n") if
                               i.startswith(b"ATOM  ") and int(i[22:26]) == mut.loc and
                               AA_SYN[mut.origAA] == i[16:20].decode('utf-8').strip()]

        if len(relavent_confidence) < 1:
            # try to find in sequence assumes reference length of 10
            if not mut.ref_seqs:
                warn_if(self._v, 2, f"Failed to find residue {mut.origAA} in {mut.loc} -- no ref seqs")
                return -1
            sequence = self._obtain_seq(resp)
            ref = next(iter(mut.ref_seqs.values()))
            index = sequence.find(ref)
            if index == -1:
                warn_if(self._v, 2, f"Faild to find residue {mut.origAA} in {mut.loc} -- can't find reference in sequence")
                return -1
            # if reference sequence is not 10 need to change "5"
            relavent_confidence = [float(i[61:66]) for i in resp.content.split(b"\n") if
                                   i.startswith(b"ATOM  ") and int(i[22:26]) == index+6]

        return relavent_confidence[0]

    def alpha_seq(self, prot):
        resp = req.get(ALPHAFOLD_PDB_URL.format(prot.Uid), timeout=TIMEOUT)
        if not resp.ok:
            return ''
        return self._obtain_seq(resp)

    def _obtain_seq(self, bytefile):
        """
        obtaipns protein sequence from bytefile response
        """
        seqs = [row[17:].split(b' ')[2:15] for row in bytefile.content.split(b"\n") if row.startswith(b"SEQRES ")]
        seqs = [AA_SYN_REV[item.decode('utf-8')] for sublist in seqs for item in sublist if item != b'']
        return "".join(seqs)

    @staticmethod
    def remove_whitespaces(text):
        """
        removes all whitespaces from string
        :param str:
        :return:
        """
        return re.sub('\s', '', text)

    def download_pdb(self, id, path):
        """
        downloads specific pdb file and save to path as id.pdb
        :param id: str id of pdb file/cs/usr/gal_passi/dina-lab/sourc
        :param path: str directory to save the file to
        :return:
        """
        self.pdpl.retrieve_pdb_file(id, pdir=path, file_format="pdb" )



