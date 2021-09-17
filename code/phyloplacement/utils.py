"""
Functions for general purposes
"""

import os


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
                         extension: str = None) -> str:
    """
    Get default path to outputfile
    """
    basename, ext = os.path.splitext(input_path)
    if extension is None:
        extension = ext
    if tag is None:
        tag = ''
    return f'{basename}{tag}{extension}'