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
from utils import print_if, adaptive_chunksize, afm_iterator, warn_if, summary_df
from definitions import *
import glob
from math import ceil
import numpy as np


def all_proteins():
    """
    :return: iterable of all Protein objects in DB
    """
    for path in glob.glob(os.path.join(PROTEIN_PATH, '*')):
        prot_name = os.path.basename(path)
        yield Protein(ref_name=prot_name)


def all_mutations():
    for prot in all_proteins():
        for mutation in prot.generate_mutations():
            yield mutation


def create_parser():
    parser = argparse.ArgumentParser(
        description="Pipeline interface for protein-level analysis of missense variants in familial data"
    )

    parser.add_argument(
        "--action",
        type=str,
        choices=["init-DB", "update-DB", "delete-DB", "score-ESM", "score-EVE", "score-AFM", "rank-DS", "to-csv"],
        default="init-DB",
        help="init-DB: initialize project database according to the supplied csv file\n"
             "score-[model]: calculate [model] scores for all missense variants in DB\n"
             "rank-DS calculates an aggregated score of all three models by normalizing and averaging available scores"
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
        "--recalc",
        type=int,
        default=0,
        choices=[0, 1],
        help="recalculate scores for all mutations including those with existing score \n"
             "1 - recalc, 0 - calculate only missing scores",
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
        default=0,
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
        if gene in REMOVED_PROTEINS:
            continue
        mut_desc = row[args.variant_col]
        dna = {'chr': row[args.chromosome_col], 'start': row[args.dna_start_col], 'end': row[args.dna_end_col],
               'ref_na': row[args.wt_col], 'alt_na': row[args.alt_col]}
        try:
            protein = Protein(ref_name=gene, verbose_level=args.verbose) if not created_protein else created_protein
            created_protein = protein
            protein.add_mut(mut_desc, dna)
        except TimeoutError:
            print_if(args.verbose, VERBOSE['thread_warnings'], f"skipped {gene} due to timeout")
            os.remove(rf'DB\test_p\{gene}')
            os.remove(rf'DB\test_m\{gene}_{mut_desc}.txt')
            skipped.append(idx)
    return skipped


def to_csv(include_type=False, outpath=''):
    df = summary_df(include_status=include_type)
    for mutation in all_mutations():
        df.loc[len(df)+1] = mutation.scores_to_csv(include_status=include_type)
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
    :param model: EVE | ESM | AFM
    :return:
    """
    assert model in AVAILABLE_MODELS, 'model should be one of EVE | ESM | AFM'
    score = None if model == 'ESM' else NO_SCORE
    score_type = NO_TYPE
    for mutation in all_mutations():
        mutation.update_score(model, score, eve_type=score_type, esm_type=score_type)
    print_if(args.verbose, VERBOSE['program_progress'], f"{model} scores set to default")


def calc_mutations_afm_scores(args, analyzer, chunk=None, iter_desc='', use_alias=False):
    """
    calculates Alpha Missense scores for all the mutations of a specific protein
    :param recalc: bool re-calculate scored for mutations with available scores
    :param use_alias: bool should reviewed uid aliases be searched
    :param args:
    :param analyzer: ProteinAnalyzer object
    :param chunk: optional df chunk to search mutation in
    :param mutations: Iterator of Mutations
    :return:
    """
    successful = 0
    tasks = [mut for mut in all_mutations() if not mut.has_afm]
    if not tasks:
        return successful
    n_muts = len(tasks)
    if chunk is not None:
        uid_index = set(chunk['uniprot_id'].unique())
    for mutation in tqdm.tqdm(tasks, desc=iter_desc, total=n_muts):
        print_if(args.verbose, VERBOSE['thread_progress'], f"Calculating AlphaMissense scores for {mutation.long_name}")
        score = analyzer.score_mutation_afm(mutation, chunk=chunk, uid_index=uid_index, use_alias=use_alias)
        successful += score is not None
        mutation.update_score('AFM', score)
        if successful == n_muts:
            return successful
    return successful


def calc_mutations_eve_scores(args, analyzer, recalc=False, iter_desc='', impute=True):
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
        score, prediction = analyzer.score_mutation_evemodel(protein=mutation.protein, mutation=mutation)
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
            score, eve_type = analyzer.score_mutation_eve_impute(mut=mutation, offset=args.offset, gz=True)
            if score != -1:
                successful += 1
                mutation.update_score('EVE', score, eve_type=eve_type)
    print_if(args.verbose, VERBOSE['program_progress'], f"done, scored {successful} of {total_tasks} mutations")
    return successful


def calc_mutations_esm_scores(args, analyzer, recalc=False, iter_desc=''):
    """
    :param args:
    :param analyzer: ProteinAnalyzer object
    :param recalc: bool re-calculate scored for mutations with available scores
    :param iter_desc:
    :return:
    """
    successful = 0
    if recalc:
        tasks = list(all_mutations())
    else:
        tasks = [mut for mut in all_mutations() if not mut.has_esm]
    total_tasks = len(tasks)
    for mutation in tqdm.tqdm(tasks, desc=iter_desc, total=len(tasks)):
        print_if(args.verbose, VERBOSE['thread_progress'], f"Calculating ESM1b scores for {mutation.long_name}")
        score, score_type = analyzer.score_mutation_esm(mut=mutation, offset=args.offset)
        if score is not None:
            successful += 1
            mutation.update_score('ESM', score, esm_type=score_type)
    print_if(args.verbose, VERBOSE['program_progress'], f"done, scored {successful} of {total_tasks} mutations")
    return successful


def calc_mutations_dsrank(n_scores_thr):
    """

    :param n_scores_thr: int number of score under which ds-rank will be considered None
    :return:
    """
    df = to_csv()
    #  MIN-MAX SCORES NORMALIZATION
    df[EVE_COL] = (df[EVE_COL] - df[EVE_COL].min()) / (df[EVE_COL].max() - df[EVE_COL].min())
    df[AFM_COL] = (df[AFM_COL] - df[AFM_COL].min()) / (df[AFM_COL].max() - df[AFM_COL].min())
    df[ESM_COL] = (df[ESM_COL] - df[ESM_COL].min()) / (df[ESM_COL].max() - df[ESM_COL].min())
    df[ESM_COL] = 1 - df[ESM_COL]
    df['n_scores'] = df[[EVE_COL, ESM_COL, AFM_COL]].count(axis=1, numeric_only=False)
    df[DS_COL] = (df[EVE_COL].fillna(0) + df[ESM_COL].fillna(0) + df[AFM_COL].fillna(0)) / df['n_scores']
    for _, row in df.iterrows():
        protein = Protein(ref_name=row[PROT_COL])
        mutation = Mutation(f"p.{row[MUT_COL]}", protein)
        score = None if row['n_scores'] < n_scores_thr else float(row[DS_COL])
        mutation.update_score(model='DS', score=score)
    # change all values under thr to nan
    df.loc[df['n_scores'] < n_scores_thr, DS_COL] = np.nan
    df.drop('n_scores', axis=1, inplace=True)
    return df


def main(args):
    workers = args.workers if args.workers else cpu_count()
    analyzer = Analyze.ProteinAnalyzer()
    for action in args.action:
        if action == 'init-DB':
            print_if(args.verbose, VERBOSE['program_progress'], f"running on {workers} CPUs")
            assert os.path.exists(args.data_path), f"could not locate data csv file at {args.data_path}"
            target = partial(create_new_records, args)
            skipped = build_db(args, target=target, workers=workers)
            return skipped
        if action == 'score-AFM':
            chunksize = adaptive_chunksize(AFM_ROWSIZE, LOW_MEM_RAM_USAGE)
            n_muts, iter_num, total_iter, total_scores =len(glob.glob(pjoin(MUTATION_PATH, '*.txt'))), 1, \
                                                        ceil(AFM_ROWS / chunksize), 0
            if args.recalc:
                erase_mutations_scores('AFM')

            print_if(args.verbose, VERBOSE['program_progress'], f"Calculating AlphaMissense scores...")
            for chunk in afm_iterator(int(chunksize), usecols=['uniprot_id', 'protein_variant', 'am_pathogenicity']):
                total_scores += calc_mutations_afm_scores(args, analyzer, chunk, f'iter {iter_num} of {total_iter} ',
                                                     use_alias=False)
                iter_num += 1
            # Second run with aliases for mutations without scores
            print_if(args.verbose, VERBOSE['program_progress'], f"Expending search for unresolved mutations...")
            iter_num = 1
            for chunk in afm_iterator(int(chunksize), usecols=['uniprot_id', 'protein_variant', 'am_pathogenicity']):
                total_scores += calc_mutations_afm_scores(args, analyzer, chunk, f'iter {iter_num} of {total_iter}',
                                                     use_alias=True)
                iter_num += 1
            print_if(args.verbose, VERBOSE['program_progress'], f"done, scored {total_scores} of {n_muts} mutations")
        if action == 'score-EVE':
            print_if(args.verbose, VERBOSE['program_progress'], f"Calculating EVEmodel scores...")
            calc_mutations_eve_scores(args, analyzer, recalc=args.recalc, impute=args.use_cpt)
        if action == 'score-ESM':
            print_if(args.verbose, VERBOSE['program_progress'], f"Calculating ESM2b scores...")
            calc_mutations_esm_scores(args, analyzer, recalc=args.recalc)
        if action == 'rank-DS':
            print_if(args.verbose, VERBOSE['program_progress'], f"Calculating Dataset Rank values...")
            return calc_mutations_dsrank(args.ds_thr)
        if action == 'to-csv':
            print_if(args.verbose, VERBOSE['program_progress'], f"extracting scored to csv...")
            return to_csv(args.include_type, args.out_path)



if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    main(args)
