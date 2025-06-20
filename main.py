import Analyze
from Mutation import Mutation
from Protein import Protein
import pandas as pd
import os
import argparse
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool as Pool
import tqdm
from functools import partial
from utils import print_if, adaptive_chunksize, afm_iterator, warn_if, summary_df, \
    esm3_setup, esm_setup, esm_seq_logits
from definitions import *
import glob
from math import ceil
import numpy as np
import shutil


def protein_names():
    protein_names = [os.path.basename(path) for path in glob.glob(os.path.join(PROTEIN_PATH, '*'))]
    return protein_names


def all_proteins():
    """
    :return: iterable of all Protein objects in DB
    """
    for path in glob.glob(os.path.join(PROTEIN_PATH, '*')):
        prot_name = os.path.basename(path)
        if prot_name in REMOVED_PROTEINS:
            continue
        yield Protein(ref_name=prot_name)


def all_mutations():
    for prot in all_proteins():
        for mutation in prot.generate_mutations():
            yield mutation


def mutations_in_csv(data, args):
    for _, row in data.iterrows():
        if row[args.protein_col] not in REMOVED_PROTEINS:
            yield Mutation(f"p.{row[args.variant_col]}", row[args.protein_col])
        else:
            pass


def create_parser():
    parser = argparse.ArgumentParser(
        description="Pipeline interface for protein-level analysis of missense variants in familial data"
    )

    parser.add_argument(
        "--action",
        type=str,
        choices=["init-DB", "update-DB", "delete-DB", "score-ESM", "score-ESM3",  "score-EVE", "score-AFM", "rank-DS",
                 "rank-family-DS", "to-csv", 'score-INTERFACE', 'score-PIONEER'],
        default="init-DB",
        help="init-DB: initialize project database according to the supplied csv file\n"
             "score-[model]: calculate [model] scores for all missense variants in DB\n"
             "rank-DS calculates an aggregated score of all three models by normalizing and averaging available scores"
             "rank-family-DS same as DS per family instead of the entire dataset"
             "use --ds-thr to set threshold for number of models required to calc scores\n"
             "to-csv outputs a summary file",
        nargs="+",
    )

    parser.add_argument(
        "--data-path",
        type=str,
        default='data.csv',
        help="path to csv file containing the data see example.csv for expected format",
    )

    parser.add_argument(
        "--out-path",
        type=str,
        default='summary.csv',
        help="path to output csv file",
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=0,
        help="number of CPUs to run on, if not given will use maximum CPUs available",
    )

    parser.add_argument(
        "--family-summary",
        type=int,
        default=0,
        choices=[0, 1],
        help="output summary csv file per family \n"
             "1 - result summary using family files, 0 - result summary for all mutations",
    )

    parser.add_argument(
        "--protein-col",
        type=str,
        default="Protein",
        help="column name holding the protein reference name in csv file i.e. PT53",
    )

    parser.add_argument(
        "--variant-col",
        type=str,
        default="Variant",
        help="column name holding the missense variant description in csv file i.e. D17Y | p.D17Y",
    )

    parser.add_argument(
        "--patient-col",
        type=str,
        default="Patient",
        help="column name holding the patient i.d",
    )

    parser.add_argument(
        "--family-col",
        type=str,
        default="Family",
        help="column name holding the family i.d",
    )

    parser.add_argument(
        "--chromosome-col",
        type=str,
        default="Chr",
        help="column name holding the chromosome in which the gene is found",
    )

    parser.add_argument(
        "--dna-start-col",
        type=str,
        default="Start",
        help="column name holding the variant dna start index",
    )

    parser.add_argument(
        "--dna-end-col",
        type=str,
        default="End",
        help="column name holding the variant dna end index",
    )

    parser.add_argument(
        "--wt-col",
        type=str,
        default="Ref",
        help="column name holding the variant dna wt (reference) nucleic acid",
    )

    parser.add_argument(
        "--alt-col",
        type=str,
        default="Alt",
        help="column name holding the variant dna altered nucleic acid",
    )

    parser.add_argument(
        "--ncbi-col",
        type=str,
        default="NCBI",
        help="column name holding the variant ncbi transcript (protein) id i.e. NM_[]",
    )

    parser.add_argument(
        "--dsrank-col",
        type=str,
        default=DS_COL,
        help="column name where ds-rank will be stored"
    )

    parser.add_argument(
        "--use-ncbi",
        type=int,
        default=1,
        choices=[0, 1],
        help="bool should NCBI_ID be used to obtain protein isoforms\n"
             "1 - use NCBI_ID, 0 - do not use NCBI_ID",
    )

    parser.add_argument(
        "--recalc",
        type=int,
        default=0,
        choices=[0, 1],
        help="recalculate scores for all mutations including those with existing score \n"
             "1 - recalc, 0 - calculate only missing scores",
    )

    parser.add_argument(
        "--unreviewed",
        type=int,
        default=0,
        choices=[0, 1],
        help="try to find scores using unreviewed entries for proteins with no scores\n"
             "currently only employed in AFM"
             "1 - use unreviewed entries, 0 - use only reviewed entries",
    )
    parser.add_argument(
        "--optimized",
        type=int,
        default=0,
        choices=[0, 1],
        help="Use optimized version for the method if available \n"
             "1 - uses optimized version of the algorithms if available, 0 - uses conventional version",
    )
    parser.add_argument(
        "--use-cpt", "--cpt",
        type=int,
        default=1,
        choices=[0, 1],
        help="use CPT to impute score when EVEmodel is not available \n"
             "1 - use cpt, 0 - only original EVEmodel",
    )

    parser.add_argument(
        "--offset",
        type=int,
        default=1,
        help="offset for mutation index default is 0",
    )

    parser.add_argument(
        "--include-type",
        type=int,
        default=0,
        choices=[0, 1],
        help="whether summary csv  should include esm and eve score types\n"
             "1 - include, 0 - omit, default",
    )

    parser.add_argument(
        "--ds-thr",
        type=int,
        default=2,
        choices=list(range(1, len(AVAILABLE_MODELS))),
        help="Number of models scores required to calculate DSRank \n"
             "default is 2",
    )

    parser.add_argument(
        "-v", "--verbose",
        type=int,
        choices=[0, 1, 2, 3],
        default=1,
        help="level of description (2-3 not recommended with workers > 1) : \n"
             "0 - Only major errors will be printed no progress information will be printed\n"
             "1 - default, program progress information and warnings will be printed\n"
             "2 - high level process and warnings i.e. protein being created\n"
             "3 - low level process information and raw errors messages i.e. detailed steps of protein creation, "
             "mutation creation etc...\n "
    )

    parser.add_argument(
        "--token",
        type=str,
        default="",
        help="Hugging-Face token required to use esm3",
    )

    parser.add_argument(
        "--esm-inference",
        type=int,
        default=0,
        choices=[0, 1],
        help="should esm scores be infered from model or used precomputed scores \n"
             "1 - inference, 0 - precomputed",
    )

    parser.add_argument(
        "--scoring-method",
        type=str,
        choices=["wt_marginals", "mutant_marginals", "masked_marginals"],
        default="masked_marginals",
        help="scoring method for esm model inference",
    )

    parser.add_argument("--task_id", type=int, required=False, help="Slurm array task ID")
    return parser


def create_new_records(args, *rowiters):
    """
    creates protein and mutation object from csv row
    :param args: user arguments
    :param rowiter: DataFrame rows iterator
    :return: bool success
    """
    idxs = [t[0] for t in rowiters]
    rows = [t[1] for t in rowiters]
    rowiter = zip(idxs, rows)
    skipped = []
    created_protein = None
    for idx, row in rowiter:
        gene = row[args.protein_col]
        ncbi_id = row[args.ncbi_col] if args.use_ncbi else ''
        if gene in REMOVED_PROTEINS:
            continue
        mut_desc = row[args.variant_col]
        dna = {'chr': row[args.chromosome_col], 'start': row[args.dna_start_col], 'end': row[args.dna_end_col],
               'ref_na': row[args.wt_col], 'alt_na': row[args.alt_col]}
        try:
            protein = Protein(ref_name=gene, ncbi=ncbi_id, verbose_level=args.verbose) if not created_protein else created_protein
            created_protein = protein
            #  same protein may have multiple ncbi ids
            if ncbi_id not in created_protein.isoforms:
                created_protein = created_protein.add_isoform(ncbi_id, 'ncbi')
            protein.add_mut(mut_desc, dna)
        except TimeoutError:
            print_if(args.verbose, VERBOSE['thread_warnings'], f"skipped {gene} due to timeout")
            os.remove(rf'DB\test_p\{gene}')
            os.remove(rf'DB\test_m\{gene}_{mut_desc}.txt')
            skipped.append(idx)
    return skipped


def to_csv(include_type=False, outpath='', mutations=None, family_name=''):
    df = summary_df(include_status=include_type)
    mutations_iter = all_mutations() if mutations is None else mutations
    for mutation in mutations_iter:
        df.loc[len(df)+1] = mutation.scores_to_csv(include_status=include_type, family_name=family_name)
    df.replace(0, np.nan, inplace=True)
    df.replace(-1, np.nan, inplace=True)
    if outpath:
        df.to_csv(outpath)
    return df


def build_db(args, target, workers):
    skipped = []
    df = pd.read_csv(args.data_path)
    values_iter = lambda value, data: data[data[args.protein_col] == value].iterrows()

    # tasks are per-protein iterators to prevent race conditions
    # may happen with proteins having multiple mutations
    unique_rows = df[~df[args.protein_col].duplicated(keep=False)]
    unique_proteins = unique_rows[args.protein_col].unique()
    tasks_unique = [values_iter(value, unique_rows) for value in unique_proteins]

    repeating_rows = df[df[args.protein_col].duplicated(keep=False)]
    repeating_proteins = repeating_rows[args.protein_col].unique()
    tasks_repeating = [values_iter(value, repeating_rows) for value in repeating_proteins]

    tasks = tasks_repeating + tasks_unique
    print_if(args.verbose, VERBOSE['program_progress'], f"Building protein database...")
    with Pool(workers) as p:
        for status in p.starmap(target, tasks):
            if status:  # if failed will return row index
                skipped += status

    if args.verbose >= VERBOSE['program_progress']:
        if not skipped:
            print_if(args.verbose, VERBOSE['program_progress'], f"done, successfully created all records")
        else:
            print_if(args.verbose, VERBOSE['program_progress'], f"done, skipped records in indexes:\n"
                                                                f"{', '.join(str(x) for x in skipped)}")

    return skipped


def erase_mutations_scores(model):
    """
    reset model scores to default
    :param model: EVE | ESM | AFM | INTERFACE
    :return:
    """
    assert model in AVAILABLE_MODELS, 'model should be one of EVE | ESM | AFM | INTERFACE'
    score = None if model in ['ESM', 'INTERFACE'] else NO_SCORE
    score_type = NO_TYPE
    for mutation in all_mutations():
        mutation.update_score(model, score, eve_type=score_type, esm_type=score_type, afm_type=score_type)
    print_if(args.verbose, VERBOSE['program_progress'], f"{model} scores set to default")


def calc_mutations_afm_scores(args, analyzer, chunk=None, iter_desc='', use_alias=True, use_unreviewed=False):
    """
    calculates Alpha Missense scores for all the mutations of a specific protein
    :param use_alias: bool if True will use uids aliases in search
    :param use_unreviewed: bool if True will use unreviewed uids aliases in search
    :param args:
    :param analyzer: ProteinAnalyzer object
    :param chunk: optional df chunk to search mutation in
    :return:
    """
    successful = 0
    tasks = [mut for mut in all_mutations() if not mut.has_afm]
    if not tasks:
        return successful
    n_muts = len(tasks)
    if chunk is not None:
        uid_index = list(chunk['uniprot_id'].unique())
        edges = uid_index[0], uid_index[-1]
    else:
        uid_index, edges = None, None
    for mutation in tqdm.tqdm(tasks, desc=iter_desc, total=n_muts):
        print_if(args.verbose, VERBOSE['thread_progress'], f"Calculating AlphaMissense scores for {mutation.long_name}")
        score, afm_type = analyzer.score_mutation_afm(mutation, chunk=chunk, uid_index=uid_index,
                                                      explicitly_access=edges, use_alias=use_alias,
                                                      use_unreviewed_uids=use_unreviewed)
        successful += score is not None
        mutation.update_score('AFM', score=score, afm_type=afm_type)
        if successful == n_muts:
            return successful
    return successful


def calc_mutations_eve_scores(args, analyzer, recalc=False, iter_desc='', impute=True, optimized=0):
    """
    :param args:
    :param analyzer: ProteinAnalyzer object
    :param recalc: bool re-calculate scored for mutations with available scores
    :param iter_desc:
    :param impute: bool should CPT be used to calculate scores when EVEmodel scores are not available
    :return:
    """
    # first iteration calculated using original EVEmodel
    successful = 0

    if recalc:
        tasks = list(all_mutations())
    else:
        tasks = [mut for mut in all_mutations() if not mut.has_eve]
    total_tasks = len(tasks)
    for mutation in tqdm.tqdm(tasks, desc=iter_desc, total=len(tasks)):
        print_if(args.verbose, VERBOSE['thread_progress'], f"Calculating EVEModel scores for {mutation.long_name}")
        score, prediction = analyzer.score_mutation_evemodel(protein=mutation.protein, mutation=mutation, optimized=optimized)
        if score != -1:
            successful += 1
            mutation.update_score('EVE', score, eve_type='EVEmodel')

    # second iteration calculated using original CPT for unresolved variants
    if impute:
        if recalc:
            tasks = [mut for mut in all_mutations() if not mut.eve_type == 'EVEmodel']
        else:
            tasks = [mut for mut in all_mutations() if not mut.has_eve]
        for mutation in tqdm.tqdm(tasks, desc=iter_desc, total=len(tasks)):
            print_if(args.verbose, VERBOSE['thread_progress'], f"Calculating CPT scores for {mutation.long_name}")
            score, eve_type = analyzer.score_mutation_eve_impute(mut=mutation, offset=args.offset, gz=True, optimized=optimized)
            if score != -1:
                successful += 1
                mutation.update_score('EVE', score, eve_type=eve_type)
    print_if(args.verbose, VERBOSE['program_progress'], f"done, scored {successful} of {total_tasks} mutations")
    return successful


def calc_mutations_esm_scores(args, analyzer, recalc=False, iter_desc='', inference=False, method='masked_marginals', optimized=0):
    """
    :param args:
    :param analyzer: ProteinAnalyzer object
    :param recalc: bool re-calculate scored for mutations with available scores
    :param iter_desc:
    :param inference: bool use precomputed scores or infer from model
    :param method: scoring method str: wt_marginals | mutant_marginals | masked_marginals
    :return:
    """
    successful = 0
    if recalc:
        tasks = list(all_mutations())
    else:
        tasks = [mut for mut in all_mutations() if not mut.has_esm]
    total_tasks = len(tasks)
    if inference:
        print_if(args.verbose, VERBOSE['program_progress'], f"Loading models...")
        model, alphabet = esm_setup(model_name=ESM1B_MODEL)
        print_if(args.verbose, VERBOSE['program_progress'], f"done")
    for mutation in tqdm.tqdm(tasks, desc=iter_desc, total=len(tasks)):
        print_if(args.verbose, VERBOSE['thread_progress'], f"Calculating ESM1b scores for {mutation.long_name}")
        if inference:
            score, score_type = analyzer.score_mutation_esm_inference(model=model, alphabet=alphabet, mut=mutation,
                                                                      method=method, offset=args.offset, log='infer\t')
        else:
            score, score_type = analyzer.score_mutation_esm1b_precomputed(mut=mutation, offset=args.offset, optimized=optimized)
        if score is not None:
            successful += 1
            mutation.update_score('ESM', score, esm_type=score_type)
    print_if(args.verbose, VERBOSE['program_progress'], f"done, scored {successful} of {total_tasks} mutations")
    return successful


def calc_mutations_dsrank(args, data_path='', out_path='', update_mut=True):
    """

    :param n_scores_thr: int number of score under which ds-rank will be considered None
    :return:
    """
    df = pd.read_csv(args.data_path) if not data_path else pd.read_csv(data_path)
    #  MIN-MAX SCORES NORMALIZATION
    df[EVE_COL] = (df[EVE_COL] - df[EVE_COL].min()) / (df[EVE_COL].max() - df[EVE_COL].min())
    df[AFM_COL] = (df[AFM_COL] - df[AFM_COL].min()) / (df[AFM_COL].max() - df[AFM_COL].min())
    df[ESM_COL] = (df[ESM_COL] - df[ESM_COL].min()) / (df[ESM_COL].max() - df[ESM_COL].min())
    df[ESM_COL] = 1 - df[ESM_COL]
    df['n_scores'] = df[[EVE_COL, ESM_COL, AFM_COL]].count(axis=1, numeric_only=False)
    df[args.dsrank_col] = (df[EVE_COL].fillna(0) + df[ESM_COL].fillna(0) + df[AFM_COL].fillna(0)) / df['n_scores']
    if update_mut:
        for _, row in df.iterrows():
            protein = Protein(ref_name=row[PROT_COL])
            mutation = Mutation(f"p.{row[MUT_COL]}", protein)
            score = None if row['n_scores'] < args.ds_thr else float(row[args.dsrank_col])
            mutation.update_score(model='DS', score=score)
    # change all values under thr to nan
    df.loc[df['n_scores'] < args.ds_thr, DS_COL] = np.nan
    df.drop('n_scores', axis=1, inplace=True)
    if out_path:
        df.to_csv(out_path)
    elif args.out_path:
        df.to_csv(args.out_path)
    return df


def is_slurm_cluster():
    return any(var.startswith("SLURM_") for var in os.environ) or shutil.which("sinfo") is not None


def calc_mutations_esm3_scores_protein(args, prot_name, recalc):
    analyzer = Analyze.ProteinAnalyzer()
    token = args.token if args.token else HUGGINGFACE_TOKEN
    model, tokenizer = esm3_setup(model_name=ESM3_MODEL, token=token)
    prot = Protein(ref_name=prot_name)
    if recalc:
        tasks = [mut for mut in prot.generate_mutations()]
    else:
        tasks = [mut for mut in prot.generate_mutations() if not mut.has_esm3]
    total_tasks, successful = len(tasks), 0
    for mutation in tqdm.tqdm(tasks, total=len(tasks)):
        score, log = analyzer.score_mutation_esm3(model=model, tokenizer=tokenizer, mut=mutation)
        if score is not None:
            successful += 1
            mutation.update_score('ESM3', float(score), esm_type=log)
    return successful


def calc_mutations_esm3_scores(args, analyzer, recalc):
    print_if(args.verbose, VERBOSE['program_progress'], f"Model inference performed on {DEVICE}")
    token = args.token if args.token else HUGGINGFACE_TOKEN
    assert token, 'must provide HuggingFace read permission token to use ESM3'
    print('loading model...')
    model, tokenizer = esm3_setup(model_name=ESM3_MODEL, token=token)
    if recalc:
        tasks = list(all_mutations())
    else:
        tasks = [mut for mut in all_mutations() if not mut.has_esm3]
    total_tasks, successful = len(tasks), 0
    for mutation in tqdm.tqdm(tasks, total=len(tasks)):
        print_if(args.verbose, VERBOSE['thread_progress'], f"Calculating ESM3 scores for {mutation.long_name}")
        score, log = analyzer.score_mutation_esm3(model=model, tokenizer=tokenizer, mut=mutation)
        if score is not None:
            successful += 1
            mutation.update_score('ESM3', float(score), esm_type=log)


def calc_mutations_interface_scores(args, analyzer, use_alias=True, use_unreviewed=False):
    """
    calculates Interface scores for all mutations
    :param use_alias: bool if True will use uids aliases in search
    :param use_unreviewed: bool if True will use unreviewed uids aliases in search
    :param args:
    :param analyzer: ProteinAnalyzer object
    :return:
    """
    successful = 0
    tasks = [mut for mut in all_mutations() if not mut.has_interface]
    if not tasks:
        return successful
    n_muts = len(tasks)
    for mutation in tqdm.tqdm(tasks, total=n_muts):
        print_if(args.verbose, VERBOSE['thread_progress'], f"Calculating Interface scores for {mutation.long_name}")
        score = analyzer.score_interface(mutation=mutation, use_alias=use_alias,
                                                         use_unreviewed_uids = use_unreviewed)
        successful += score is not None
        mutation.update_score('INTERFACE', score=score)
        if successful == n_muts:
            return successful
    return successful


def main(args):
    workers = args.workers if args.workers else cpu_count()
    analyzer = Analyze.ProteinAnalyzer()
    n_muts = len(glob.glob(pjoin(MUTATION_PATH, '*.txt')))
    for action in args.action:
        if action == 'init-DB':
            print_if(args.verbose, VERBOSE['program_progress'], f"running on {workers} CPUs")
            assert os.path.exists(args.data_path), f"could not locate data csv file at {args.data_path}"
            target = partial(create_new_records, args)
            skipped = build_db(args, target=target, workers=workers)
            return skipped
        if action == 'score-AFM':
            #  chunk search is now depreciated
            #chunksize = adaptive_chunksize(AFM_ROWSIZE, LOW_MEM_RAM_USAGE)
            #n_muts, iter_num, total_iter, total_scores =len(glob.glob(pjoin(MUTATION_PATH, '*.txt'))), 1, ceil(AFM_ROWS / chunksize), 0
            if args.recalc:
                erase_mutations_scores('AFM')
            print_if(args.verbose, VERBOSE['program_progress'], f"Calculating AlphaMissense scores...")
            #for chunk in afm_iterator(int(chunksize), usecols=['uniprot_id', 'protein_variant', 'am_pathogenicity']):
            #    total_scores += calc_mutations_afm_scores(args, analyzer, chunk, f'iter {iter_num} of {total_iter} ',
            #                                         use_alias=False)
            #    iter_num += 1
            total_scores = calc_mutations_afm_scores(args, analyzer, use_alias=True)
            if args.unreviewed:
                print_if(args.verbose, VERBOSE['program_progress'], f"using unreviewed entries for unresolved variants")
                total_scores += calc_mutations_afm_scores(args, analyzer, use_alias=True, use_unreviewed=True)
            print_if(args.verbose, VERBOSE['program_progress'], f"done, scored {total_scores} mutations")
        if action == 'score-EVE':
            print_if(args.verbose, VERBOSE['program_progress'], f"Calculating EVEmodel scores...")
            calc_mutations_eve_scores(args, analyzer, recalc=args.recalc, impute=args.use_cpt, optimized=args.optimized)
        if (action == 'score-ESM') or (action =='score-ESM1b'):
            print_if(args.verbose, VERBOSE['program_progress'], f"Calculating ESM-1b scores...")
            calc_mutations_esm_scores(args, analyzer, recalc=args.recalc, inference=args.esm_inference,
                                      method=args.scoring_method, optimized=args.optimized)
        if action =='score-ESM3':
            print_if(args.verbose, VERBOSE['program_progress'], f"Calculating ESM-3 scores...")
            if is_slurm_cluster():
                if args.optimized == 0:
                    calc_mutations_esm3_scores(args, analyzer, recalc=args.recalc)
                else:
                    proteins = protein_names()
                    if args.task_id <= len(proteins):
                        prot_name = proteins[args.task_id - 1]
                        print(f"Processing {args.task_id}: {prot_name}")
                        calc_mutations_esm3_scores_protein(args=args, prot_name=prot_name, recalc=args.recalc)
                    else:
                        print(f"Task ID {args.task_id} exceeds available protein count ({len(proteins)}). Skipping.")
            else:
                calc_mutations_esm3_scores(args, analyzer, recalc=args.recalc)
        if action in ['score-INTERFACE', 'score-PIONEER']:
            if args.recalc:
                erase_mutations_scores('INTERFACE')
            print_if(args.verbose, VERBOSE['program_progress'], f"Calculating interface scores...")
            total_scores = calc_mutations_interface_scores(args, analyzer, use_alias=True)
            if args.unreviewed:
                print_if(args.verbose, VERBOSE['program_progress'], f"using unreviewed entries for unresolved variants")
                total_scores += calc_mutations_interface_scores(args, analyzer, use_alias=True, use_unreviewed=True)
            print_if(args.verbose, VERBOSE['program_progress'], f"done, scored {total_scores} of {n_muts} mutations")
        if action == 'rank-DS':
            print_if(args.verbose, VERBOSE['program_progress'], f"Calculating Dataset Rank values...")
            calc_mutations_dsrank(args)
            print_if(args.verbose, VERBOSE['program_progress'], f"done")
        if action == 'rank-family-DS':
            print_if(args.verbose, VERBOSE['program_progress'], f"Calculating Family Datasets Rank values...")
            family_files_path = glob.glob(pjoin(FAMILY_SCORES_PATH, '*.csv'))
            for path in tqdm.tqdm(family_files_path, total=len(family_files_path)):
                family_name = os.path.basename(path)[:-4]
                outpth = pjoin(FAMILY_SCORES_PATH, f"famrank_{family_name}.csv")
                calc_mutations_dsrank(args, data_path=path, out_path=outpth, update_mut=False)
            print_if(args.verbose, VERBOSE['program_progress'], f"done")
        if action == 'to-csv':
            print_if(args.verbose, VERBOSE['program_progress'], f"extracting scored to csv...")
            if args.family_summary:
                family_files_path = glob.glob(pjoin(FAMILY_DATA_PATH, '*.csv'))
                for path in tqdm.tqdm(family_files_path, total=len(family_files_path)):
                    name = os.path.basename(path)[:-4]
                    family_data = pd.read_csv(path)
                    mutations = mutations_in_csv(family_data, args)
                    outpth = pjoin(FAMILY_SCORES_PATH, f"famrank_{name}.csv")
                    to_csv(args.include_type, outpath=outpth, mutations=mutations, family_name=family_name)
            else:
                to_csv(args.include_type, args.out_path)
            print_if(args.verbose, VERBOSE['program_progress'], f"done")


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    main(args)
