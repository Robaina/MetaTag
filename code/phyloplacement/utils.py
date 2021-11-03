"""
Functions and classes for general purposes
"""

import os
import shutil
import random
import string
import subprocess
import tarfile
import pickle
from functools import partial
from multiprocessing import Pool


def handle_exceptions(foo):
    """
    Decorator to handle possible exceptions in
    given function (foo)
    """
    def inner_foo(*args, **kwargs):
        try:
            foo(*args, **kwargs)
        except Exception as e:
            print(f'{foo.__name__} failed with exception: {e}')
    return inner_foo
class TemporaryFilePath:
    """
    Custom context manager to create a temporary file
    which is removed when exiting context manager
    """
    def __init__(self,
                 work_dir: str = None,
                 extension: str = None,
                 create_file: bool = False):
        self.work_dir = work_dir or ''
        self.extension = extension or ''
        self.create_file = create_file

    def __enter__(self):
        temp_id = ''.join(
        random.choice(string.ascii_lowercase) for i in range(10)
        )
        self.file_path = os.path.join(
            self.work_dir, f'temp_{temp_id}{self.extension}'
            )
        if self.create_file:
            os.mkdir(self.file_path)
        return self.file_path

    def __exit__(self, exc_type, exc_val, exc_tb):
        if os.path.exists(self.file_path):
            os.remove(self.file_path)
class TemporaryDirectoryPath:
    """
    Custom context manager to create a temporary directory
    which is removed when exiting context manager
    """
    def __init__(self, work_dir: str = None):
        self.work_dir = work_dir or ''

    def __enter__(self):
        temp_id = ''.join(
        random.choice(string.ascii_lowercase) for i in range(10)
        )
        self.dir_path = os.path.join(
            self.work_dir, f'temp_{temp_id}/'
            )
        os.mkdir(self.dir_path)
        return self.dir_path

    def __exit__(self, exc_type, exc_val, exc_tb):
        if os.path.exists(self.dir_path):
            shutil.rmtree(self.dir_path)

def createTemporaryFilePath(work_dir: str = None, extension: str = None):
    """
    Converted into custom context manager
    """
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
    TODO: return absolute path
    """
    basename = os.path.basename(input_path)
    dirname = os.path.abspath(os.path.dirname(input_path))
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
        stdout = subprocess.DEVNULL
    else:
        stdout = None
    output = subprocess.run(
        command_str, shell=True,
        cwd=work_dir, capture_output=return_output,
        stdout=stdout
        )
    return output

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

def listTarDir(tar_dir: str) -> list:
    """
    List files within tar or tar.gz directory
    """
    with tarfile.open(tar_dir, 'r') as tar_obj:
        files = tar_obj.getnames()
    return files
