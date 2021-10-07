"""
Functions for general purposes
"""

import os
import pickle


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

def terminalExecute(command_str: str, suppress_output=False) -> None:
    """
    Execute given command in terminal through Python
    """
    if suppress_output:
        suppress_code = '>/dev/null 2>&1'
        command_str = f'{command_str} {suppress_code}'
    os.system(command_str)
    
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

def countRecords(fasta_file: str) -> None:
    cmd_str = f'grep -c ">" {fasta_file}'
    terminalExecute(cmd_str)

# def sliceFasta(input_file, output_file, N):
#     n = 0
#     records = SeqIO.parse(input_file, 'fasta')
#     sliced_records = []
#     for record in records:
#         if n < N:
#             sliced_records.append(record)
#         else:
#             break
#         n += 1
#     SeqIO.write(sliced_records, output_file, 'fasta')