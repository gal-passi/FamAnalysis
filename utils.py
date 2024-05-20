#import multiprocessing, geneffect, firm
import os
import glob
import warnings

import Family
import Patient
import Protein
import Mutation
import Pair

ALPHAFOLD_PDB_URL = "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v1.pdb"

def print_if(verbose, thr, text):
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

def make_fasta(path, name,seq):
    full_path = os.path.join(path, f"{name}.fasta")
    with open(full_path, "w+") as file:
        file.write(f">{name}\n")
        file.write(f"{seq}\n")


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
        cur_patients = patients[0:i] + patients[i+1:len(patients)]
        sets = [set(p.mutations()) for p in cur_patients]
        final_set |= set.intersection(*sets)
    return final_set


def recurring_family_mutations(families = ['HR1', 'HR3', 'HR4', 'HR5', 'HR6', 'HR7', 'HR8', 'HR9', 'HR10', 'HR11', 'HR12']):
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
                data[f'bert_{i+1}'].append(score)
            data['bert_sum'].append(sum(list(details['bertScore'])))
            data['bert_max'].append(max(list(details['bertScore'])))
        else:
            for i in range(1,6):
                data[f'bert_{i}'].append(-1)
            data['bert_sum'].append(-1)
            data['bert_max'].append(-1)
        data['mentioned'].append(str(Analyze.ProteinAnalyzer.search_litriture(prot.name)))


def find_mutations_with_pdbs(protein, log = "log.txt", found = 'found.txt'):
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


def all_proteins():
    """
    :return: iterable of all Protein objects in DB
    """
    for path in glob.glob(os.path.join(PROTEIN_DB, '*')):
        prot_name = os.path.basename(path)
        yield Protein(ref_name=prot_name)


def all_mutations():
    for path in glob.glob(os.path.join(MUTATION_DB, '*')):
        mut_name = os.path.basename(path)
        if mut_name == 'desktop.ini':
            continue
        yield utils.create_mutation(mut_name)


#TODO there is a problem with firm setup might not be compatiable with newer version of Bio. maybe need to update firm
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