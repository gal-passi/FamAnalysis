import pandas as pd
from utils import afm_iterator, adaptive_chunksize, ugzip
from utils import SafeDownloader
from definitions import *
import os
from os.path import join as pjoin
import json
import zipfile
import gzip
import ast
import shutil
import tqdm
from math import ceil
import numpy as np


def create_directories():
    print("building directory...")
    os.mkdir(DB_PATH)
    for subdir in DB_CHILDREN:
        os.mkdir(pjoin(DB_PATH, subdir))
    for subdir in FAMILY_SUBDIRS:
        os.mkdir(pjoin(FAMILY_PATH, subdir))
    print("done")


def create_afm_index():
    index = {}
    chuksize = adaptive_chunksize(AFM_UID_ROWSIZE, DEFAULT_RAM_USAGE)
    print('building Alpha Missense index...')
    for chunk in afm_iterator(chuksize, usecols=['uniprot_id']):
        chunk.drop_duplicates(inplace=True, keep='first')
        ind = chunk['uniprot_id'].to_dict()
        ind.update(index)  # keep first index
        index = ind
    inverse_index = {v: k for k, v in index.items()}
    with open(AFM_DIRECTORY_PATH, "w") as file:
        file.write(json.dumps(inverse_index))
    idxs = list(inverse_index.values())
    idxs.sort()
    ranges = {int(idxs[i]): int(idxs[i + 1]) for i in range(len(idxs) - 1)}
    ranges[int(idxs[-1])] = None
    with open(AFM_RANGES_PATH, "w") as file:
        file.write(json.dumps(ranges))
    print('done')
    return inverse_index, ranges


def split_afm(index, ranges, chunksize=10**6, delte_main_file=True):
    os.mkdir(pjoin(AFM_PATH, AFM_PROTEINS))
    truncated_record, truncated_uid = None, None
    print('splitting Alpha Missense main file...')
    df_iter =  pd.read_csv(AFM_DATA_PATH, sep='\t', chunksize=chunksize, header=AFM_HEADER,
                           usecols=['uniprot_id', 'protein_variant', 'am_pathogenicity'])
    for df in tqdm.tqdm(df_iter, '', ceil(AFM_ROWS / chunksize)):
        uids = df['uniprot_id'].unique()
        for uid in uids:
            record = df[df['uniprot_id'] == uid]
            idx_start = index[uid]
            idx_end = ranges[str(idx_start)]
            if idx_end is None:  # last record
                record.to_csv(fr'DB/AFM/protein_files/{uid}.csv')
                continue
            expected_length = idx_end - idx_start
            if uid == truncated_uid:
                record = pd.concat([truncated_record, record])
                assert len(record) == expected_length, f'concat error {uid}: expected:{expected_length} record: {len(record)}'
                truncated_record, truncated_uid = None, None
            elif len(record) < expected_length:
                assert truncated_uid is None, f'new truncated {uid} before resolving previous truncated {uid}'
                truncated_record = record
                truncated_uid = uid
                continue
            record.to_csv(pjoin(AFM_PROTEINS_PATH, f'{uid}.csv.gz'), compression='gzip')
    if delte_main_file:
        os.remove(AFM_DATA)


def create_esm_index():
    df = pd.read_csv(pjoin(ESM_DATA_PATH, 'contents_u_df.csv'))
    index = {row['gene']: row['id'] for _, row in df.iterrows()}
    with open(ESM_INDEX_PATH, "w") as file:
        file.write(json.dumps(index))

def process_pioneer():
    pioneer_sequences, pioneer_seq_predictions = {}, {}
    pioneer_df = pd.read_csv(PIONEER_DATA_PATH, delimiter='\t')
    for _, row in tqdm.tqdm(pioneer_df.iterrows(), total=len(pioneer_df)):
        uid_1, uid_2, source, inter_1, inter_2, seq_1, seq_2 = tuple(row)
        inter_1, inter_2 = ast.literal_eval(inter_1), ast.literal_eval(inter_2)
        for uid, seq, pred in zip([uid_1, uid_2], [seq_1, seq_2], [inter_1, inter_2]):
            if uid not in pioneer_sequences:
                pioneer_sequences[uid] = {seq}
            else:
                pioneer_sequences[uid].add(seq)
            if seq not in pioneer_seq_predictions:
                pioneer_seq_predictions[seq] = np.zeros(len(seq), dtype=int)
            #  this will only retain last prediction source
            if pred:
                pioneer_seq_predictions[seq][np.array(pred) - 1] = SOURCE_TO_NUM[source]
    #  serialize data for json
    pioneer_sequences = {k: list(v) for k, v in pioneer_sequences.items()}
    pioneer_seq_predictions = {k: v.tolist() for k, v in pioneer_seq_predictions.items()}
    os.chdir(PIONEER_PATH)
    with open('pioneer_sequences.json', 'w') as file:
        file.write(json.dumps(pioneer_sequences, indent=4))
    with open('pioneer_predictions.json', 'w') as file:
        file.write(json.dumps(pioneer_seq_predictions, indent=4))
    os.chdir(ROOT_DIR)
    print('done')


def download_afm():
    os.chdir(AFM_PATH)
    downloader = SafeDownloader([AFM_PUBLIC_DATA], [AFM_DATA + '.gz'])
    print('Downloading AlphaMissense source data...')
    downloader.download()
    print('Unzipping AlphaMissense data...')
    chunksize = adaptive_chunksize(AFM_ROWSIZE, DEFAULT_RAM_USAGE)
    ugzip(path=AFM_DATA + '.gz', outfile=AFM_DATA, chunksize=chunksize)
    os.chdir(ROOT_DIR)
    print('done')


def download_eve(remove_zip=True):
    os.chdir(EVE_PATH)
    if EVE_PARTIAL_DOWNLOAD:
        urls = [EVE_SINGLE_PROTEIN.format(name) for name in EVE_PARTIAL_DOWNLOAD]
        names = [name + '.zip' for name in EVE_PARTIAL_DOWNLOAD]
        downloader = SafeDownloader(urls, names)
    else:
        downloader = SafeDownloader([EVE_PUBLIC_DATA], [EVE_DATA + '.zip'])
    print('Downloading EVEModel source data...')
    downloader.download()
    print('Unzipping EVEModel data...')
    if EVE_PARTIAL_DOWNLOAD:
        for name in names:
            with zipfile.ZipFile(name + '.zip', "r") as zip_ref:
                zip_ref.extractall(name+'.csv')
    else:
        with zipfile.ZipFile(EVE_DATA + '.zip', "r") as zip_ref:
            zip_ref.extractall(EVE_DATA)
        if remove_zip:
            os.remove(EVE_DATA + '.zip')
        os.chdir(EVE_DATA_PATH)
        for file in os.listdir('.'):
            if os.path.basename(file) != EVE_VARIANTS:
                shutil.rmtree(file)
    os.chdir(ROOT_DIR)
    shutil.move(pjoin(EVE_SETUP_INDEX_DIR, 'eve_index.txt'), EVE_INDEX_PATH_2)
    shutil.move(pjoin(EVE_SETUP_INDEX_DIR, 'eve_reverse_index.txt'), EVE_INVERSE_INDEX)
    shutil.rmtree(EVE_SETUP_INDEX_DIR)
    print('done')


def download_cpt(remove_zip=True):
    os.chdir(CPT_PATH)
    filenames = [CPT_EVE_DATA_NAME + '.zip', CPT_IMPUTE_DATA_1_NAME + '.zip',
                                 CPT_IMPUTE_DATA_2_NAME + '.zip']
    downloader = SafeDownloader([CPT_EVE_DATA, CPT_IMPUTE_DATA_1, CPT_IMPUTE_DATA_2], filenames)
    print('Downloading CPT source data...')
    downloader.download()
    print('Unzipping CPT data...')
    for name in filenames:
        with zipfile.ZipFile(name, "r") as zip_ref:
            zip_ref.extractall(name[:-4])
    # MERGE EVE IMPUTE FOLDERS
    shutil.copytree(CPT_EVE_IMPUTE_PATH_2, CPT_EVE_IMPUTE_PATH_1, dirs_exist_ok=True)
    os.rename(CPT_IMPUTE_DATA_1_NAME, CPT_IMPUTE_DATA_NAME)
    shutil.rmtree(CPT_IMPUTE_DATA_2_NAME)
    os.chdir(CPT_PATH)
    if remove_zip:
        for file in filenames:
            os.remove(file)
    os.chdir(ROOT_DIR)
    print('done')


def download_esm(remove_zip=True):
    os.chdir(ESM_PATH)
    downloader = SafeDownloader([ESM_VARIANTS_DATA], [ESM_DATA + '.zip'])
    print('Downloading ESM-1b variants source data...')
    downloader.download()
    print('Unzipping ESM-1b data...')
    with zipfile.ZipFile(ESM_DATA + '.zip', "r") as zip_ref:
        zip_ref.extractall(ESM_DATA)
    if remove_zip:
        os.remove(ESM_DATA + '.zip')
    create_esm_index()
    os.chdir(ROOT_DIR)
    print('done')

def download_pioneer():
    os.chdir(PIONEER_PATH)
    downloader = SafeDownloader([PIONEER_PUBLIC_DATA], [PIONEER_RAW_DATA])
    print('Downloading PIONEER data...')
    downloader.download()
    os.chdir(ROOT_DIR)
    print('done')

def download_all_data():
    download_afm()
    download_eve()
    download_cpt()
    download_esm()
    download_pioneer()


def main():
    create_directories()
    download_afm()
    create_afm_index()
    download_esm()
    create_esm_index()
    download_eve()
    download_cpt()
    download_pioneer()
    process_pioneer()
    print("Setup completed successfully")


if __name__ == '__main__':
    main()