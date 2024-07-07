import glob

from utils import progress_bar, afm_iterator, adaptive_chunksize, ugzip
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

def download_afm():
    os.chdir(AFM_PATH)
    downloader = SafeDownloader([AFM_PUBLIC_DATA], [AFM_DATA])
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
    print('done')

def download_cpt(remove_zip=True):
    os.chdir(CPT_PATH)
    filenames = [CPT_EVE_DATA_NAME + '.zip', CPT_IMPUTE_DATA_1_NAME + '.zip',
                                 CPT_IMPUTE_DATA_2_NAME + '.zip']
    downloader = SafeDownloader([CPT_EVE_DATA, CPT_IMPUTE_DATA_1, CPT_IMPUTE_DATA_2], filenames)
    print('Downloading CPT source data...')
    downloader.download()
    print('Unzipping EVEModel data...')
    for name in filenames:
        with zipfile.ZipFile(name, "r") as zip_ref:
            zip_ref.extractall(name[:-4])
    # MERGE EVE IMPUTE FOLDERS
    shutil.copytree(CPT_EVE_IMPUTE_PATH_2, CPT_EVE_IMPUTE_PATH_1, dirs_exist_ok=True)
    os.rename(CPT_IMPUTE_DATA_1_NAME, CPT_IMPUTE_DATA_NAME)
    shutil.rmtree(CPT_IMPUTE_DATA_2_NAME)
    if remove_zip:
        for file in filenames:
            os.remove(file)
    os.chdir(ROOT_DIR)
    print('done')

def download_all_data():
    download_afm()
    download_eve()
    download_cpt()


if __name__ == '__main__':
    #create_directories()
    #download_data()
    #create_afm_index()
    os.chdir(CPT_PATH)
    shutil.rmtree('CPT1_score_no_EVE_set_2')
    os.chdir(ROOT_DIR)

