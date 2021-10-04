"""
Tools to preprocess sequence databases

1. Parse and clean file paths
2. Parse and reformat sequence labels and assert sequences are in right format
"""

import os
import sys
import re
import shutil
import pyfastx

from .utils import saveToPickleFile, setDefaultOutputPath

upper_lower_digits = re.compile('[^a-zA-Z0-9]')
upper_case_letters = re.compile('[^A-Z]')


def reIndexFASTA(input_fasta: str, output_dir: str = None):
    """
    Change record ids for numbers and store then in a dictionary
    Useful when converting fasta alignments to phylip format 
    (phylip ids are less than 10 char long)
    """
    if output_dir is None:
        output_dir = os.path.dirname(input_fasta)
    
    id_dict = {}
    fasta_file = setDefaultOutputPath(input_fasta, 
                                      tag='_short_ids',
                                      only_filename=True)
    dict_file = setDefaultOutputPath(input_fasta,
                                     tag='_id_dict',
                                     extension='.pickle',
                                     only_filename=True)
    output_fasta = f'{os.path.join(output_dir, fasta_file)}'
    output_dict = f'{os.path.join(output_dir, dict_file)}'
    
    fa = pyfastx.Fasta(input_fasta)
    id_dict = dict(zip(range(len(fa)), fa.keys()))
    with open(output_fasta, 'w') as fp:
        for record_id in fa.keys():
            seq = fa[record_id]
            fp.write(seq.raw.replace(seq.description, str(seq.id - 1)))

    saveToPickleFile(id_dict, output_dict)

def cleanFilePath(file_name: str) -> None:
    """
    Asserts file path has correct format
    """
    fdir = os.path.dirname(file_name)
    fname, ext = os.path.splitext(os.path.basename(file_name))
    clean_fname = upper_lower_digits.sub(
        '_', fname).replace('__', '_').strip('_')
    return os.path.join(fdir, f'{clean_fname}{ext}')

def assertPeptideSequenceIsLegit(record_seq) -> str:
    """
    Assert peptide sequence is legit
    """
    return upper_case_letters.sub('', record_seq.upper())

def modifyRecordID(record_name: str, record_number: int,
                   id_type: int,
                   file_name: str = None) -> str:
    """
    id_type: '1: number', '2: description', '3: filename_number'
    """
    if id_type == 1:
        return record_number
    elif id_type == 2:
        return record_name
    else:
        return f'{file_name}_{record_number}'

def reformatFASTAfile(fasta_file: str, id_type: int,
                      output_file: str = None) -> None:
    """
    Reformat FASTA file ok.
    """
    dirname = os.path.dirname(fasta_file)
    basename = os.path.basename(fasta_file)
    fname, ext = os.path.splitext(basename)

    if output_file is None:
        output_file = os.path.join(dirname, f'{fname}_modified{ext}')
    else:
        output_file = os.path.abspath(output_file)

    fasta = pyfastx.Fasta(fasta_file, build_index=False, full_name=True)
    with open(output_file, 'w') as outfile:
        n = 0
        for record_name, record_seq in fasta:
            n += 1
            record_id = modifyRecordID(record_name,
                                       record_number=n,
                                       id_type=id_type,
                                       file_name=fname)
            record_seq = assertPeptideSequenceIsLegit(record_seq)

            outfile.write(f'>{record_id}\n{record_seq}\n')

def pipe_line(fasta_path: str, id_type: int,
              output_file = None) -> None:
    """
    Pipeline!!
    """
    def is_legit_path(fasta_path, legit_fasta_path):
        return fasta_path == legit_fasta_path
            
    clean_fasta_path = cleanFilePath(fasta_path)

    if not is_legit_path(fasta_path, clean_fasta_path):
        shutil.move(fasta_path, clean_fasta_path)
    
    reformatFASTAfile(fasta_file=clean_fasta_path,
                      id_type=id_type, output_file=output_file)


if __name__ == '__main__':

    """
    python3 modify_fasta.py path/to/fasta [path/to/dir/] id_type

    id_type:
      1: Only numbers
      2: Long description
      3: File name + number
    """

    fasta_path = sys.argv[1]
    id_type = int(sys.argv[2])
    fasta_dir = os.path.dirname(fasta_path)

    # type_answer = int(input(
    #     ('Should I change labels for file name followed '
    #     'by a numbers (1) or leave as is (2) or just '
    #     'numbers (3)?')))

    path_is_directory = os.path.isdir(fasta_path)

    if path_is_directory:
        file_names = os.listdir(fasta_path)
        
        for file_name in file_names:
            fasta_path = os.path.join(fasta_dir, file_name)
            pipe_line(fasta_path, id_type=id_type)
    else:
        pipe_line(fasta_path, id_type=id_type)

    print('I am done!')