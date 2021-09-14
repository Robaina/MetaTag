"""
Functions to preprocess sequence data
"""

from .utils import terminalExecute


def removeDuplicates():
    """
    Remove sequences in fasta files with repeated labels
    or sequences.

    NOTES:

    script: joinseqs.py
    """
    pass

def reformatFileName():
    """
    Checks if file name format is legal
    and corrects otherwise


    NOTES:

    script: clean.py

    1. File names, must contain only upper/lower case letters and digits and '_',
    replace anything else (such as a space) by '_'

    2. Sequences must be only composed of uppercase letters A-Z
    """
    pass

def refomatSequencesAndLabels():
    """
    Checks for inconsistencies in sequence data and
    reformat sequence labels

    NOTES:

    script: clean.py

    1. File names, must contain only upper/lower case letters and digits and '_',
    replace anything else (such as a space) by '_'

    2. Sequences must be only composed of uppercase letters A-Z
    """
    pass

def reformatProdigalLabels():
    """
    Reformat prodigal labels and file names

    NOTES:

    script: reformatabel.py
    """
    pass

def translateDNA():
    """
    Translate DNA sequences with prodigal

    NOTES:

    script: loopparallel.py
    """
    pass

def filterDataByTargetProtein(hmm_model: str):
    """
    Generate protein-specific database by filtering
    sequence database to only contain sequences 
    corresponing to protein of interest

    NOTES:

    script: pfamgenomesparallel.py
    """
    pass