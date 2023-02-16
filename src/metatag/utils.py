#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Functions and classes for general purposes
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import pickle
import random
import shutil
import string
import subprocess
import tarfile
from functools import partial
from multiprocessing import Pool


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


def handle_exceptions(foo):
    """
    Decorator to handle possible exceptions in
    given function (foo)
    """

    def inner_foo(*args, **kwargs):
        try:
            foo(*args, **kwargs)
        except Exception as e:
            print(f"{foo.__name__} failed with exception: {e}")

    return inner_foo


class TemporaryFilePath:
    """
    Custom context manager to create a temporary file
    which is removed when exiting context manager
    """

    def __init__(
        self, work_dir: str = None, extension: str = None, create_file: bool = False
    ):
        self.work_dir = work_dir or ""
        self.extension = extension or ""
        self.create_file = create_file

    def __enter__(self):
        temp_id = "".join(random.choice(string.ascii_lowercase) for i in range(10))
        self.file_path = os.path.join(self.work_dir, f"temp_{temp_id}{self.extension}")
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
        self.work_dir = work_dir or ""

    def __enter__(self):
        temp_id = "".join(random.choice(string.ascii_lowercase) for i in range(10))
        self.dir_path = os.path.join(self.work_dir, f"temp_{temp_id}/")
        os.mkdir(self.dir_path)
        return self.dir_path

    def __exit__(self, exc_type, exc_val, exc_tb):
        if os.path.exists(self.dir_path):
            shutil.rmtree(self.dir_path)


def create_temporary_file_path(work_dir: str = None, extension: str = None):
    """
    Converted into custom context manager
    """
    if work_dir is None:
        work_dir = ""
    if extension is None:
        extension = ""
    temp_id = "".join(random.choice(string.ascii_lowercase) for i in range(10))
    return os.path.join(work_dir, f"temp_{temp_id}{extension}")


def delete_temporary_files(dir_path: str) -> None:
    """
    Remove files from directory
    """
    for fname in os.listdir(dir_path):
        os.remove(os.path.join(dir_path, fname))


def set_default_output_path(
    input_path: str,
    tag: str = None,
    extension: str = None,
    only_filename: bool = False,
    only_dirname: bool = False,
) -> str:
    """
    Get default path to output file or directory
    """
    basename = os.path.basename(input_path)
    dirname = os.path.abspath(os.path.dirname(input_path))
    fname, ext = os.path.splitext(basename)
    if extension is None:
        extension = ext
    if tag is None:
        tag = ""
    default_file = f"{fname}{tag}{extension}"
    if only_filename:
        return default_file
    if only_dirname:
        return os.path.abspath(dirname)
    else:
        return os.path.abspath(os.path.join(dirname, default_file))


def save_to_pickle_file(python_object, path_to_file="object.pkl"):
    """
    Save python object to pickle file
    """
    out_file = open(path_to_file, "wb")
    pickle.dump(python_object, out_file)
    out_file.close()


def read_from_pickle_file(path_to_file="object.pkl"):
    """
    Load python object from pickle file.
    Returns python object.
    """
    in_file = open(path_to_file, "rb")
    python_object = pickle.load(in_file)
    in_file.close()
    return python_object


def terminal_execute(
    command_str: str,
    suppress_shell_output=False,
    work_dir: str = None,
    return_output=False,
) -> subprocess.STDOUT:
    """
    Execute given command in terminal through Python
    """
    if suppress_shell_output:
        stdout = subprocess.DEVNULL
    else:
        stdout = None
    output = subprocess.run(
        command_str,
        shell=True,
        cwd=work_dir,
        capture_output=return_output,
        stdout=stdout,
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


def full_path_list_dir(dir: str) -> list:
    """
    Return full path of files in provided directory
    """
    return [os.path.join(dir, file) for file in os.listdir(dir)]


def extract_tar_file(tar_file: str, dest_dir: str = None) -> None:
    """
    Extract tar or tar.gz files to dest_dir
    """
    if dest_dir is None:
        dest_dir = "."
    if tar_file.endswith("tar.gz"):
        tar = tarfile.open(tar_file, "r:gz")
        tar.extractall(path=dest_dir)
        tar.close()
    elif tar_file.endswith("tar"):
        tar = tarfile.open(tar_file, "r:")
        tar.extractall(path=dest_dir)
        tar.close()
    else:
        raise ValueError("Input is not a tar file")


def list_tar_dir(tar_dir: str) -> list:
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
    def from_pickle_paths(cls, dict_paths: list[str]) -> DictMerger:
        """
        Initialize class from list of paths to dictionaries (pickle)
        """
        dict_list = [
            cls.read_from_pickle_file(dict_path.strip()) for dict_path in dict_paths
        ]
        return cls(dicts=dict_list)

    @staticmethod
    def read_from_pickle_file(path_to_file="object.pkl"):
        """
        Load python object from pickle file.
        Returns python object.
        """
        in_file = open(path_to_file, "rb")
        python_object = pickle.load(in_file)
        in_file.close()
        return python_object

    def merge(
        self, dict_prefixes: list[str] = None, save_pickle_path: str = None
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
