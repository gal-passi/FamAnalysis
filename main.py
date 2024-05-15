import Analyze
import main
import utils
from Protein import Protein
import Mutation
import Connections
from Analyze import ProteinAnalyzer
import re
import pandas as pd
import requests as req
import urllib.request
import glob
import os
import json
import pickle
import subprocess
from subprocess import check_output
import warnings


CONTACT = "gal.passi@mail.huji.ac.il"
HEADERS = {'User-Agent': 'Python {}'.format(CONTACT)}
UNIPORT_URL = "http://www.uniprot.org/uniprot/"
PROTEIN_DB = 'DB/proteins'
MUTATION_DB = 'DB/Mutations'

#data = {'name': [], 'variant': [], 'firmScore': [], 'eveScore': [], 'evePrediction': [], 'bert_1': [], 'bert_2': [], 'bert_3': [], 'bert_4': [], 'bert_5': [], 'mentioned': [], 'pdbs':[]}

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


if __name__ == "__main__":
    scores1 = []
    df = pd.read_csv('data_normalization.csv')
    for _, row in df.iterrows():
        prot = Protein(ref_name=row['protein'], LOA)
        mut = Mutation.Mutation(f"p.{row['variant']}", prot)
        if mut.bert_score != 0:
            scores1.append(mut.bert_score)
    print(scores1)

'''
if mut.eve_score != -1 or mut.eve_score != '-1':
    print("has eve")
    prot.mutations[mut.extended_description]['eve_impute'] = False
else:
    score = a.score_mutation_eve_impute(mut)
    print(score)
    prot.mutations[mut.extended_description]['eve_impute'] = True
    prot.mutations[mut.extended_description]['eveScore'] = score
    print("==================")

prot._update_DB(os.path.join(prot.directory, prot.MUTS), prot.mutations, mode='pickle')
        prot = mut.protein
prot.muts[mut.extended_description]['main_entery'] = c.entery_name(prot)
prot.muts[mut.extended_description]['all_entries'] = c.entery_name(prot, all_results=True)
'''

'''
c = Connections.Uniport()
    analyzer = ProteinAnalyzer()
    n_prots = len(list(glob.glob(os.path.join(PROTEIN_DB, '*'))))
    processed = 0
    for prot in all_proteins():
        processed += 1
        if not prot.mutations:
            continue
        for variant, details in prot.mutations.items():
            mut = Mutation.Mutation(variant, prot, details)
            print(f"processing {prot.name}({mut.name})")
            res = analyzer.score_mutation_evemodel(prot, mut)
            if res != (-1, -1):
                print(f"has EVE {res}")
                prot.mutations[variant]['eveScore'] = res[0]
                prot.mutations[variant]['evePrediction'] = res[1]
                prot.mutations[mut.extended_description]['eve_impute'] = False
            else:
                score = analyzer.score_mutation_eve_impute(mut, gz=True)
                print(f'CPT score {score}')
                prot.mutations[mut.extended_description]['eve_impute'] = True
                prot.mutations[mut.extended_description]['eveScore'] = score
        prot._update_DB(os.path.join(prot.directory, prot.MUTS), prot.mutations, mode='pickle')
        if processed % 30 == 0:
            print(f"{processed} / {n_prots}")
'''
'''
    a = Analyze.ProteinAnalyzer()
    c = Connections.Uniport()
    for prot in all_proteins():
        print(prot.name)
        prot._Uids['main_entery'] = c.entery_name(prot)
        prot._Uids['all_enteries'] = c.entery_name(prot, all_results=True)
        with open(os.path.join(prot.directory, prot.UIDS), "w") as file:
            file.write(json.dumps(prot._Uids))
'''

'''
#UPDATE MUTATION
df = pd.read_csv("varients_segregation_distribution.csv")
for protein in all_proteins():
    print(f"processing {protein.name}:")
    muts = protein.mutations
    temp = {mut: {'chr': 0, 'start': 0, 'end': 0, "firmScore": -1, "eveScore": -1, "bertScore": -1} for mut in muts}
    for mut in muts:
        for index, data in df.iterrows():
            if mut in data['AAChange_refGeneMT']:
                temp[mut]['chr'] = data['Chr']
                temp[mut]['start'] = int(data['Start'])
                temp[mut]['end'] = int(data['End'])

    protein._update_DB(os.path.join(protein.directory, protein.MUTS), temp, mode='pickle')
'''

'''         
#EVEMODEL

    processed = 0
    _unip = Connections.Uniport()
    analyzer = ProteinAnalyzer()
    for mut in all_mutations():
        processed += 1
        if processed % 25 == 0:
            print(f"processed {processed}...")
        if not mut.eve_score != -1:
            continue
        else:
            print(mut.long_name)
            prot = mut.protein
            desc = mut.extended_description
            res = analyzer.score_mutation_evemodel(prot, mut)
            prot.mutations[desc]['eveScore'] = res[0]
            prot.mutations[desc]['evePrediction'] = res[1]
            prot._update_DB(os.path.join(prot.directory, prot.MUTS), prot.mutations, mode='pickle')

'''


'''
#FIRM
    analyzer = ProteinAnalyzer()
    firm_classifier, geneffect_setup, thread_pool = analyzer._firm_setup()
    for prot in all_proteins():
        if not prot.mutations:
            continue
        else:
            for desc, dna_data in prot.mutations.items():
                if dna_data['firmScore'] == -1.0:
                    mut = Mutation.Mutation(desc, prot, dna_data)
                    print(f"processing {prot.name}({mut.name})")
                    try:
                        snp = geneffect_setup.variant_interpreter.process_snp(mut.chr, mut.start, mut.ref_na, mut.alt_na)
                        snp_gene_effect, = snp.cds_gene_effects
                        prot.mutations[desc]['firmScore'] = firm_classifier.predict_adjusted_proba(snp_gene_effect)
                        print("done")
                    except Exception as e:
                        print("failed")
                        print("*********************")
                        print(str(e))
                        print("*********************")
                        prot.mutations[desc]['firmlog'] = str(e)

            prot._update_DB(os.path.join(prot.directory, prot.MUTS), prot.mutations, mode='pickle')
'''

'''
#UPDATE MUTATION
df = pd.read_csv("varients_segregation_distribution.csv")
for protein in all_proteins():
    print(f"processing {protein.name}:")
    muts = protein.mutations
    temp = {mut: {'chr': 0, 'start': 0, 'end': 0, "firmScore": -1, "eveScore": -1, "bertScore": -1} for mut in muts}
    for mut in muts:
        for index, data in df.iterrows():
            if mut in data['AAChange_refGeneMT']:
                temp[mut]['chr'] = data['Chr']
                temp[mut]['start'] = int(data['Start'])
                temp[mut]['end'] = int(data['End'])

    protein._update_DB(os.path.join(protein.directory, protein.MUTS), temp, mode='pickle')
    print("done")
'''

'''
#BERT CSV TO PROTEIN
    bert_records = list(map(lambda x: os.path.basename(x), glob.glob("DB/bert/" + "*.csv")))
    bert_records = [i for i in bert_records if '_' in i]
    bert_records = [i for i in bert_records if 'iso' not in i]
    total = len(bert_records)
    done = 0
    for record in bert_records:
        prot_name = record.split("_")[0]
        mut_name = record.split("_")[1][:-4]
        df = pd.read_csv(f"DB/bert/{record}")
        try:
            data = (float(df['esm1v_t33_650M_UR90S_1']), float(df['esm1v_t33_650M_UR90S_2']), float(df['esm1v_t33_650M_UR90S_3']),
                    float(df['esm1v_t33_650M_UR90S_4']), float(df['esm1v_t33_650M_UR90S_5']))
        except KeyError:
            print("problem (bert) with:")
            print(prot_name)
            print(mut_name)
            print("-"*60)
            done += 1
            continue

        prot = Protein(ref_name=prot_name)
        try:
            if prot.mutations[mut_name]['bertScore']:
                done += 1
                continue
            prot.mutations[mut_name]['bertScore'] = data
        except KeyError:
            print("problem (mutation) with:")
            print(prot_name)
            print(mut_name)
            print("-"*60)
            done += 1
            continue
        prot._update_DB(os.path.join(prot.directory, prot.MUTS), prot.mutations, mode='pickle')
        done += 1
        if done % 25 == 0:
            print(f"finished {round(done/total, 2)}...")
'''




'''
SUSPEICOUS GENES
for prot in all_proteins():
        for mname, data in prot.mutations.items():
            escore = data['eveScore']
            fscore = data['firmScore'] if data['firmScore'] not in ['err', 'asr', 'asr_1'] else -1
            if (escore >= 0.3) and (escore < 0.7):
                inter_eve.add(f"{prot.name} - {mname} - {escore} - {data['evePrediction']}")
            if escore >= 0.7:
                severe_eve.add(f"{prot.name} - {mname} - {escore} - {data['evePrediction']}")
            if (fscore >= 0.3) and (fscore < 0.65):
                inte_firm.add(f"{prot.name} - {mname} - {fscore}")
            if (fscore >= 0.65) and (fscore < 0.8):
                severe_firm.add(f"{prot.name} - {mname} - {fscore}")
            if fscore >= 0.8:
                critical_firm.add(f"{prot.name} - {mname} - {fscore}")
    print("inter_eve:")
    for data in inter_eve:
        print(data)
    print("-"*70)
    print("severe_eve:")
    for data in severe_eve:
        print(data)
    print("-"*70)
    print("inte_firm:")
    for data in inte_firm:
        print(data)
    print("-"*70)
    print("severe_firm:")
    for data in severe_firm:
        print(data)
    print("-"*70)
    print("critical_firm:")
    for data in critical_firm:
        print(data)
    print("-"*70)
'''

'''
#create csv
    analyzer = Analyze.ProteinAnalyzer()
    n_prots = len(list(glob.glob(os.path.join(PROTEIN_DB, '*'))))
    processed = 0
    for prot in all_proteins():
        processed += 1
        if processed % 25 == 0:
            print(f"finished {processed / n_prots}...")
        if not prot.mutations:
            continue
        for desc, details in prot.mutations.items():
            mut = Mutation(desc, prot, details)
            data['name'].append(prot.name)
            data['variant'].append(desc)
            data['firmScore'].append(details['firmScore'])
            data['eveScore'].append(details['eveScore'])
            data['evePrediction'].append(details['evePrediction'])
            data['mentioned'].append(str(analyzer.search_litriture(prot.name)[1]))
            data['pdbs'].append(mut.pdbs)
            
            if (details['bertScore']) and (details['bertScore'] != -1):
                data['bert_1'].append(details['bertScore'][0])
                data['bert_2'].append(details['bertScore'][1])
                data['bert_3'].append(details['bertScore'][2])
                data['bert_4'].append(details['bertScore'][3])
                data['bert_5'].append(details['bertScore'][4])
            else:
                data['bert_1'].append(-1)
                data['bert_2'].append(-1)
                data['bert_3'].append(-1)
                data['bert_4'].append(-1)
                data['bert_5'].append(-1)
            
    print(data)
    df = pd.DataFrame(data)
    df.to_csv("summary1505.csv")
'''

'''
#FULL PIPELINE
   c = Connections.Uniport()
   analyzer = ProteinAnalyzer()
   firm_classifier, geneffect_setup, thread_pool = analyzer._firm_setup()
   n_prots = len(list(glob.glob(os.path.join(PROTEIN_DB, '*'))))
   processed = 0
   bert_flag=False
   for prot in all_proteins():
       processed += 1
       if not prot.mutations:
           continue
       for variant, details in prot.mutations.items():
           mut = Mutation.Mutation(variant, prot, details)
           print(f"processing {prot.name}({mut.name})")
           if details['firmScore'] in [-1.0, 'asr', 'err']:
               try:
                   snp = geneffect_setup.variant_interpreter.process_snp(mut.chr, mut.start, mut.ref_na, mut.alt_na)
                   snp_gene_effect, = snp.cds_gene_effects
                   prot.mutations[variant]['firmScore'] = firm_classifier.predict_adjusted_proba(snp_gene_effect)
               except Exception as e:;
                   print("firm failed")
                   prot.mutations[variant]['firmlog'] = str(e)
           if details['eveScore'] == -1:
               res = analyzer.score_mutation_evemodel(prot, mut)
               prot.mutations[variant]['eveScore'] = res[0]
               prot.mutations[variant]['evePrediction'] = res[1]

           if (not details['bertScore']) or (details['bertScore'] == -1):

               analyzer.Bert_score(prot)
               record = f"{prot.name}_{variant}"
               try:
                   df = pd.read_csv(f"DB/bert/{record}.csv")
               except FileNotFoundError:
                   print("no file created by bert")
                   continue
               try:
                   bert_data = (float(df['esm1v_t33_650M_UR90S_1']), float(df['esm1v_t33_650M_UR90S_2']), float(df['esm1v_t33_650M_UR90S_3']),
                           float(df['esm1v_t33_650M_UR90S_4']), float(df['esm1v_t33_650M_UR90S_5']))
                   prot.mutations[variant]['bertScore'] = bert_data
               except KeyError:
                   print("problem with bert (partial file created)")
                   pass
'''
'''
#INTERFACE SCORE
    flag = False
    dataset = pd.read_csv('summary2405.csv')
    problems = []
    all_muts = zip(dataset['name'].to_list(), dataset['variant'].to_list())
    n = len(list(all_muts))
    c = 0
    for pname, mname in zip(dataset['name'].to_list(), dataset['variant'].to_list()):
        c += 1
        if (pname == 'BRD2') and (mname == 'p.A105V'):
            flag = True
        if flag:
            print(f"procesing {pname}({mname})")
            try:
                prot = Protein(ref_name=pname, load_only=True)
            except NameError:
                problems.append((pname, mname))
                continue
            try:
                mut = Mutation.Mutation(mname, prot)
            except NameError:
                problems.append((pname, mname))
                continue
            temp = mut.interface
            mut._save_obj(mut._directory)
            if c % 100 == 0:
                print(f'FINISHED {round(c/n,2)}%')
    print(problems)
'''

'''
#SEND PAIRS TO ALPHFOLD:
    analyzer = ProteinAnalyzer()
    with open(r'datasets/patients_pairs.txt', 'rb') as file:
        pairs = pickle.load(file)
    for patient, data in pairs.items():
        for pair in data:
            m1,m2 = pair
            mut1, mut2 = Mutation.Mutation(m1[1], m1[0]), Mutation.Mutation(m2[1], m2[0])
            analyzer.model_pairs_alphafold(mut1, mut2)

# MODEL STRUCTURE
    analyzer = ProteinAnalyzer()
    fastas = glob.glob(os.path.join('temp_alphfold_jobs', '*'))
    n = len(fastas)
    i = 0
    for file in fastas:
        i += 1
        os.system(f"python3.7 ~dina/disk/collabFold_script.py {file}")
        print(f"finished {i}/{n}")
'''

'''
for idx, row in df.iterrows():
    prot = Protein.Protein(ref_name = row['name'])
    temp = prot.muts
    temp[row['variant']]['consensusScore'] = float(row['avg_score'])
    prot._update_DB(os.path.join(prot.directory, prot.MUTS), temp, mode='pickle')
'''

'''
# MANUAL BERT
    a = Analyze.ProteinAnalyzer()
    seq = 'SSSGTSILTGSAIQVQNIKKDQTLKARIEIPSCKDVAPVEKTIKLLPSSHVARLQIFSVEGQKAIQIKHQDEVNWIAGDIMHNLIFQMYDEGEREINITSALAEKIKVNWTPEINKEHLLQGLLPDVQVPTSVKDMRYCQVSFQDDHVSLESAFTVRPLPDEPKHLKCEMKGGKTVQMGQELQGEVVIIITDQYGNQIQAFSPSSLSSLSIAGVGLDSSNLKTTFQENTQSISVRGIKFIPGPPGNKDLCFTWREFSDFIRVQLISGPPAKLLLIDWPELKESIPVINGRDLQNPIIVQLCDQWDNPAPVQHVKISLTKASNLKLMPSNQQHKTDEKGRANLGVFSVFAPRGEHTLQVKAIYNKSIIEGPIIKLMILPDPEKPVRLNVKYDKDASFLAGGLFTDFMISVISEDDSIIKNINPARISMKMWKLSTSGNRPPANAETFSCNKIKDNDKEDGCFYFRDKVIPNKVGTYCIQFGFMMDKTNILNSEQVIVEVLPNQPVKLVPKIKPPTPAVSNVRSVASRTLVRDLHLSITDDYDNHTGIDLVGTIIATIKGSNEEDTDTPLFIGKVRTLEFPFVNGSAEIMSLVLAESSPGRDSTEYFIVFEPRLPLLSRTLEPYILPFMFYNDVKKQQQMAALTKEKDQLSQSIVMYKSLFEASQQLLNEMKCQVEEARLKEAQLRNELKIHNIDIPTTQQVPHIEALLKRKLSEQEELKKKPRRSCTLPNYTKGSGDVLGKIAHLAQIEDDRAAMVISWHLASDMDCVVTLTTDAARRIYDETQGRQQVLPLDSIYKKTLPDWKRSLPHFRNGKLYFKPIGDPVFARDLLTFPDNVEHCETVFGMLLGDTIILDNLDAANHYRKEVVKITHCPTLLTRDGDRIRSNGKFGGLQNKAPPMDKLRGMVFGAPVPKQCLILGEQIDLLQQYRSAVCKLDSVNKDLNSQLEYLRTPDMRKKKQELDEHEKNLKLIEEKLGMTPIRKCNDSLRHSPKVETTDCPVPPKRMRREATRQNRIITKTDV'
    df = pd.DataFrame({'mutation': ['R1972C']})
    path = r'/cs/labs/dina/gal_passi/breast_cancer/sources/pdbs/SMCHD1_R1972C.csv'
    df.to_csv(path)
    print(f'csv created at {path}')
    script_location = os.path.join("/cs/labs/dina/gal_passi/breast_cancer/pdbs-venv/esm", "variant-prediction/predict.py")
    with subprocess.Popen(f"python3 {script_location} "
                          f"--model-location esm1v_t33_650M_UR90S_1 esm1v_t33_650M_UR90S_2 esm1v_t33_650M_UR90S_3 esm1v_t33_650M_UR90S_4 esm1v_t33_650M_UR90S_5 "
                          f"--sequence {seq} "
                          f"--dms-input {path} "
                          f"--mutation-col mutation "
                          f"--dms-output {path} "
                          f"--offset-idx {986} "
                          f"--scoring-strategy wt-marginals", stdout=subprocess.PIPE, shell=True) as proc:
        print(proc.stdout.read())


'''

'''
afm_scores = []
afm_class = []
for _, row in data.iterrows():
    afm_prot = df[df['uniprot_id'] == row['uid']]
    if afm_prot.empty:
        afm_scores.append(-1)
        afm_class.append(-1)
        continue
    afm_variant = afm_prot[afm_prot['protein_variant'] == row['variant']]
    if afm_variant.empty:
        afm_scores.append(-1)
        afm_class.append(-1)
        continue
    afm_scores.append(float(afm_variant['am_pathogenicity']))
    afm_class.append(afm_variant['am_class'])

'''
'''
# ADD NEW PROTEINS 
    df = pd.read_csv('DB/Rare_segragating_variants_families1_12.AF.01-02.csv')
    skipped = []
    for _, row in df.iterrows():

            gene = row['Gene']
            mut_desc = row['AA']
            dna = {'chr': row['Chr'], 'start': row['Start'] , 'end' : row['End'], 'ref_na': row['Ref'], 'alt_na': row['Alt']}
            try:
                protein = Protein(ref_name=gene)
                protein.add_mut(mut_desc, dna)
                print("===========================")
            except TimeoutError:
                print(f"skipped {gene} due to timeout")
                skipped.append(gene)
                os.remove(rf'DB\proteins_rare\{gene}')
                os.remove(rf'DB\Mutations_rare\{gene}_{mut_desc}.txt')
                continue


'''