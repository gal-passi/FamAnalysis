from utils import progress_bar, afm_iterator, adaptive_chunksize, ugzip
from definitions import *
import os
from os.path import join as pjoin
import wget
import json
import pandas as pd


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

def download_data():
    print('Downloading AlphaMissense source data...')
    os.chdir(AFM_PATH)
    wget.download(url=AFM_PUBLIC_DATA, bar=progress_bar, out=AFM_PATH)
    print('Unzipping AlphaMissense data...')
    chunksize = adaptive_chunksize(AFM_ROWSIZE, DEFAULT_RAM_USAGE)
    ugzip(path=AFM_DATA + '.gz', outfile=AFM_DATA, chunksize=chunksize)
    os.chdir(ROOT_DIR)
    print('done')

if __name__ == '__main__':
    #create_directories()
    #download_data()
    create_afm_index()

