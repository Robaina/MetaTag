"""
Functions for general purposes
"""

import os
import random
import string
import subprocess
import tarfile
import pickle
from functools import partial
from multiprocessing import Pool


def handle_exceptions(foo):
    def inner_foo(*args, **kwargs):
        try:
            foo(*args, **kwargs)
        except Exception as e:
            print(f'{foo.__name__} failed with exception: {e}')
    return inner_foo

def saveToPickleFile(python_object, path_to_file='object.pkl'):
    """
    Save python object to pickle file
    """
    out_file = open(path_to_file,'wb')
    pickle.dump(python_object, out_file)
    out_file.close()
    
def readFromPickleFile(path_to_file='object.pkl'):
    """
    Load python object from pickle file.
    Returns python object.
    """
    in_file = open(path_to_file,'rb')
    python_object = pickle.load(in_file)
    return python_object

def terminalExecute(command_str: str,
                    suppress_shell_output=False,
                    work_dir: str = None,
                    return_output=False) -> subprocess.STDOUT:
    """
    Execute given command in terminal through Python
    """
    if suppress_shell_output:
        suppress_code = '>/dev/null 2>&1'
        command_str = f'{command_str} {suppress_code}'
    output = subprocess.run(
        command_str, shell=True,
        cwd=work_dir, capture_output=return_output)
    return output

def createTemporaryFilePath(work_dir: str = None, extension: str = None):
    if work_dir is None:
        work_dir = ''
    if extension is None:
        extension = ''
    temp_id = ''.join(
        random.choice(string.ascii_lowercase) for i in range(10)
        )
    return os.path.join(work_dir, f'temp_{temp_id}{extension}')
    
def deleteTemporaryFiles(dir_path: str) -> None:
    """
    Remove files from directory
    """
    for fname in os.listdir(dir_path):
        os.remove(os.path.join(dir_path, fname))

def setDefaultOutputPath(input_path: str, tag: str = None,
                         extension: str = None,
                         only_filename: bool = False,
                         only_dirname: bool = False) -> str:
    """
    Get default path to outputfile
    """
    basename = os.path.basename(input_path)
    dirname = os.path.dirname(input_path)
    fname, ext = os.path.splitext(basename)
    if extension is None:
        extension = ext
    if tag is None:
        tag = ''
    default_file = f'{fname}{tag}{extension}'
    if only_filename:
        return default_file
    if only_dirname:
        return dirname
    else:
        return os.path.join(dirname, default_file)

def parallelizeOverInputFiles(callable,
                              input_list: list,
                              n_processes: int = None,
                              **callable_kwargs) -> None: 
    """
    Parallelize callable over a set of input objects using a pool 
    of workers. Inputs in input list are passed to the first argument
    of the callable.
    Additional callable named arguments may be passed.
    """
    if n_processes is None:
        n_processes = os.cpu_count - 1
    p = Pool(processes=n_processes)
    p.map(partial(callable, **callable_kwargs), input_list)
    p.close()
    p.join()

def fullPathListDir(dir: str) -> list:
    """
    Return full path of files in provided directory
    """
    return [os.path.join(dir, file) for file in os.listdir(dir)]

def extractTarFile(tar_file: str, dest_dir: str = None) -> None:
    """
    Extract tar or tar.gz files to dest_dir
    """ 
    if dest_dir is None:
        dest_dir = '.'
    if tar_file.endswith('tar.gz'):
        tar = tarfile.open(tar_file, 'r:gz')
        tar.extractall(path=dest_dir)
        tar.close()
    elif tar_file.endswith('tar'):
        tar = tarfile.open(tar_file, 'r:')
        tar.extractall(path=dest_dir)
        tar.close()
    else:
        raise ValueError('Input is not a tar file')