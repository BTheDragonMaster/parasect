# -*- coding: utf-8 -*-

"""Helpers module for PARASECT."""

import os
import shutil
from typing import List


def clear_temp_dir(dir_path: str, keep: List[str]) -> None:
    """Clear the contents of a directory, except for specified files or directories.

    :param dir_path: Path to directory to clear.
    :type dir_path: str
    :param keep: List of file names or directory names to keep.
    :type keep: List[str]
    :raises NotADirectoryError: If the directory does not exist.
    """
    # check if directory exists
    if not os.path.isdir(dir_path):
        raise NotADirectoryError(f"directory does not exist: {dir_path}")

    for file_name in os.listdir(dir_path):

        # remove all content from the temp directory, except for files in the keep
        # list (e.g., ".gitkeep")
        if file_name in keep:
            continue

        # construct path to file or directory
        path = os.path.join(dir_path, file_name)

        # remove file or directory
        if os.path.isfile(path):
            os.remove(path)

        elif os.path.isdir(path):
            shutil.rmtree(path)
