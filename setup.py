import pandas as pd
from utils import afm_iterator, adaptive_chunksize, ugzip
from utils import SafeDownloader
from definitions import *
import os
from os.path import join as pjoin
import json
import zipfile
import shutil


def create_directories():
    print("building directory...")
    os.mkdir(DB_PATH)
    for subdir in DB_CHILDREN:
        os.mkdir(pjoin(DB_PATH, subdir))
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


def create_esm_index():
    df = pd.read_csv(pjoin(ESM_DATA_PATH, 'contents_u_df.csv'))
    index = {row['gene']: row['id'] for _, row in df.iterrows()}
    with open(ESM_INDEX_PATH, "w") as file:
        file.write(json.dumps(index))


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


def download_all_data():
    download_afm()
    download_eve()
    download_cpt()
    download_esm()


def main():
    create_directories()
    download_afm()
    create_afm_index()
    download_esm()
    create_esm_index()
    download_eve()
    download_cpt()
    print("Setup completed successfully")


if __name__ == '__main__':
    main()