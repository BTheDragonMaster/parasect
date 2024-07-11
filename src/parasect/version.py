# -*- coding: utf-8 -*-

"""Version information for :mod:`parasect`.

Run with ``python -m parasect.version``
"""

import os
from subprocess import CalledProcessError, check_output  # noqa: S404

__all__ = [
    "VERSION",
    "get_version",
    "get_git_hash",
]

VERSION = "1.0.0-dev"


def get_git_hash() -> str:
    """Get the :mod:`parasect` git hash.

    :return: The git hash.
    :rtype: str
    """
    with open(os.devnull, "w", encoding="utf-8") as devnull:
        try:
            ret = check_output(  # noqa: S603,S607
                ["git", "rev-parse", "HEAD"],
                cwd=os.path.dirname(__file__),
                stderr=devnull,
            )
        except CalledProcessError:
            return "UNHASHED"
        else:
            return ret.strip().decode("utf-8")[:8]


def get_version(with_git_hash: bool = False) -> str:
    """Get the :mod:`parasect` version string, including a git hash.

    :param with_git_hash: Include the git hash in the version string.
    :type with_git_hash: bool
    :return: The version string.
    :rtype: str
    """
    return f"{VERSION}-{get_git_hash()}" if with_git_hash else VERSION


if __name__ == "__main__":
    print(get_version(with_git_hash=True))  # noqa:T201
