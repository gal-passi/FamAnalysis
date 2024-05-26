from definitions import *
import os
from os.path import dirname, abspath
from os.path import join as pjoin

def create_directories(verbose=1):
    print("building directory...")
    os.mkdir(DB_PATH)
    for subdir in DB_CHILDREN:
        os.mkdir(pjoin(DB_PATH, subdir))
    print("done")

if __name__ == '__main__':
    create_directories()


