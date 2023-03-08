#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions and classes for general purposes
"""

from __future__ import annotations

import json
import logging
import os
import pickle
import random
import shutil
import string
import subprocess
import sys
import tarfile
from functools import partial
from multiprocessing import Pool
from pathlib import Path

logger = logging.getLogger(__name__)


class ConfigParser:
    """Handle MetaTag configuration file."""

    def __init__(self, config_file: Path) -> None:
        self._config_file = Path(config_file)
        self._config = self.get_config()

    @classmethod
    def get_default_config(cls):
        """Initialize ConfigParser with default config file."""
        return cls(cls.initialize_config_file())

    @staticmethod
    def initialize_config_file() -> Path:
        """Initialize empty config file.
        Returns:
            Path: path to generated config file.
        """
        config_file = Path(__file__).parent / "config.json"
        if not config_file.exists():
            config = {
                "input_database": "",
                "hmm_directory": "",
                "max_hmm_reference_size": "",
                "min_sequence_length": "",
                "max_sequence_length": "",
                "output_directory": "",
                "translate": "",
                "relabel": "",
                "remove_duplicates": "",
                "relabel_prefixes": "",
                "hmmsearch_args": "",
                "tree_method": "",
            }
            with open(config_file, "w", encoding="UTF-8") as f:
                json.dump(config, f, indent=4)
        return config_file

    def get_config_path(self) -> Path:
        """Show config file path."""
        return self._config_file

    def get_config(self) -> dict:
        """Load config file.
        Returns:
            dict: dict containing fields and values of config file.
        """
        with open(self._config_file, "r", encoding="UTF-8") as file:
            config = json.loads(file.read())
        return config

    def write_config(self) -> None:
        """Write config dict to file."""
        with open(self._config_file, "w", encoding="UTF-8") as f:
            json.dump(self._config, f, indent=4)

    def update_config(self, key: str, value: str) -> None:
        """Update config file
        Args:
            key (str): config file key name to be updated.
            value (str): new value.
        """
        self._config[key] = value
        self.write_config()

    def get_field(self, key: str) -> str:
        """Get field from config file.
        Args:
            key (str): key name to get the value from.
        Returns:
            str: key value.
        """
        return self._config[key]


class TemporaryFilePath:
    """
    Custom context manager to create a temporary file
    which is removed when exiting context manager
    """

    def __init__(
        self, work_dir: Path = None, extension: str = None, create_file: bool = False
    ):
        self.work_dir = Path(work_dir) or ""
        self.extension = extension or ""
        self.create_file = create_file

    def __enter__(self):
        temp_id = "".join(random.choice(string.ascii_lowercase) for i in range(10))
        temp_file_name = f"temp_{temp_id}{self.extension}"
        self.file_path = (
            self.work_dir / temp_file_name if self.work_dir else Path(temp_file_name)
        )
        if self.create_file:
            self.file_path.mkdir(parents=True, exist_ok=True)
        return self.file_path

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file_path.exists():
            self.file_path.unlink()


class TemporaryDirectoryPath:
    """
    Custom context manager to create a temporary directory
    which is removed when exiting context manager
    """

    def __init__(self, work_dir: Path = None):
        self.work_dir = Path(work_dir) or ""

    def __enter__(self):
        temp_id = "".join(random.choice(string.ascii_lowercase) for i in range(10))
        self.dir_path = (
            self.work_dir / f"temp_{temp_id}/"
            if self.work_dir
            else Path(f"temp_{temp_id}/")
        )
        self.dir_path.mkdir(parents=True, exist_ok=True)
        return self.dir_path

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.dir_path.exists():
            shutil.rmtree(self.dir_path)


def create_temporary_file_path(work_dir: Path = None, extension: str = None):
    """
    Converted into custom context manager
    """
    if work_dir is None:
        work_dir = ""
    else:
        work_dir = Path(work_dir)
    if extension is None:
        extension = ""
    temp_id = "".join(random.choice(string.ascii_lowercase) for i in range(10))
    return work_dir / f"temp_{temp_id}{extension}"


def delete_temporary_files(dir_path: Path) -> None:
    """
    Remove files from directory
    """
    dir_path = Path(dir_path)
    for fname in dir_path.iterdir():
        fname.unlink()


def set_default_output_path(
    input_path: Path,
    tag: str = None,
    extension: str = None,
    only_filename: bool = False,
    only_basename: bool = False,
    only_dirname: bool = False,
) -> Path:
    """Utility function to generate a default path to output file
    or directory based on an input file name and path.
    Args:
        input_path (Path): path to input file.
        tag (str, optional): text tag to be added to file name. Defaults to None.
        extension (str, optional): change input file extension with this one. Defaults to None.
        only_filename (bool, optional): output only default filename. Defaults to False.
        only_basename (bool, optional): output only default basename (no extension). Defaults to False.
        only_dirname (bool, optional): output only path to default output directory. Defaults to False.
    Returns:
        Path: a path or name to a default output file.
    """
    input_path = Path(input_path)
    dirname = input_path.parent
    fname, ext = input_path.stem, input_path.suffix
    if extension is None:
        extension = ext
    if tag is None:
        tag = ""
    default_file = f"{fname}{tag}{extension}"
    if only_basename:
        return Path(fname)
    if only_filename:
        return Path(default_file)
    if only_dirname:
        return dirname
    else:
        return (dirname / default_file).absolute()


def save_to_pickle_file(python_object: object, path_to_file: Path = "object.pkl"):
    """
    Save python object to pickle file
    """
    outfile = open(path_to_file, "wb")
    pickle.dump(python_object, outfile)
    outfile.close()


def read_from_pickle_file(path_to_file: Path = "object.pkl"):
    """
    Load python object from pickle file.
    Returns python object.
    """
    infile = open(path_to_file, "rb")
    python_object = pickle.load(infile)
    infile.close()
    return python_object


def terminal_execute(
    command_str: str,
    suppress_shell_output=False,
    work_dir: Path = None,
    return_output=False,
) -> subprocess.STDOUT:
    """Execute given command in terminal through Python.
    Args:
        command_str (str): terminal command to be executed.
        suppress_shell_output (bool, optional): suppress shell output. Defaults to False.
        work_dir (Path, optional): change working directory. Defaults to None.
        return_output (bool, optional): whether to return execution output. Defaults to False.
    Returns:
        subprocess.STDOUT: subprocess output.
    """
    if suppress_shell_output:
        command_str = f"{command_str} >/dev/null 2>&1"
    else:
        command_str = command_str

    output = subprocess.run(
        command_str, shell=True, cwd=work_dir, capture_output=return_output, stdout=None
    )
    return output


def parallelize_over_input_files(
    callable, input_list: list, n_processes: int = None, **callable_kwargs
) -> None:
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


def extract_tar_file(tar_file: Path, dest_dir: Path = None) -> None:
    """
    Extract tar or tar.gz files to dest_dir
    """
    tar_file = Path(tar_file)
    if dest_dir is None:
        dest_dir = "."
    if tar_file.as_posix().endswith("tar.gz"):
        tar = tarfile.open(tar_file, "r:gz")
        tar.extractall(path=dest_dir)
        tar.close()
    elif tar_file.as_posix().endswith("tar"):
        tar = tarfile.open(tar_file, "r:")
        tar.extractall(path=dest_dir)
        tar.close()
    else:
        logger.error("Input is not a tar file")
        sys.exit(1)


def list_tar_dir(tar_dir: Path) -> list:
    """
    List files within tar or tar.gz directory
    """
    with tarfile.open(tar_dir, "r") as tar_obj:
        files = tar_obj.getnames()
    return files


def easy_pattern_matching(
    text: str, left_pattern: str, right_pattern: str = None
) -> str:
    """
    Just straightforward string searchs between two patterns
    """
    idx1 = text.find(left_pattern)
    left_subtext = text[idx1 + len(left_pattern) :]
    if right_pattern is not None:
        idx2 = left_subtext.find(right_pattern)
        matched_text = left_subtext[:idx2]
    else:
        matched_text = left_subtext
    return matched_text


class DictMerger:
    def __init__(self, dicts: list[dict]) -> None:
        """
        Toos to merge python dictionaries into a single one
        @param
        dicts: list of dictionaries to be merged
        """
        self._dict_list = dicts

    @classmethod
    def from_pickle_paths(cls, dict_paths: list[Path]) -> DictMerger:
        """
        Initialize class from list of paths to dictionaries (pickle)
        """
        dict_list = [
            cls.read_from_pickle_file(dict_path.as_posix().strip())
            for dict_path in dict_paths
        ]
        return cls(dicts=dict_list)

    @staticmethod
    def read_from_pickle_file(path_to_file: Path = "object.pkl"):
        """
        Load python object from pickle file.
        Returns python object.
        """
        in_file = open(path_to_file, "rb")
        python_object = pickle.load(in_file)
        in_file.close()
        return python_object

    def merge(
        self, dict_prefixes: list[str] = None, save_pickle_path: Path = None
    ) -> dict:
        """
        Merge dictionaries
        @params
        dict_prefixes: list of strings containing prefixes to be added to
                values in each dict (optional)
        """
        if dict_prefixes is None:
            dict_prefixes = ["" for _ in range(len(self._dict_list))]
        else:
            dict_prefixes = dict_prefixes

        merged_dict = {
            key: prefix + value
            for prefix, dict_object in zip(dict_prefixes, self._dict_list)
            for (key, value) in dict_object.items()
        }

        if save_pickle_path is not None:
            save_to_pickle_file(merged_dict, save_pickle_path)
        return merged_dict
