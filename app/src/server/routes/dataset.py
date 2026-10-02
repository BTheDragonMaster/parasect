# -*- coding: utf-8 -*-

"""Routes for downloading the full reference database and describing its schema."""

from __future__ import annotations

import functools
import hashlib
import os
import sqlite3

from flask import Blueprint, Response, send_file

from parasect.version import get_version

from .common import ResponseData, Status
from .constants import DB_PATH

blueprint_dataset = Blueprint("dataset", __name__)


def _download_name() -> str:
    """Name the file after the release, so a downloaded copy says which data it holds."""
    return f"parasect-{get_version()}.db"


@functools.lru_cache(maxsize=4)
def _sha256(path: str, mtime: float) -> str:
    """Checksum of the database, cached per file version (mtime is part of the key).

    :param path: Path to the database.
    :param mtime: Modification time of the database, only used as cache key.
    :return: Hex digest.
    """
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_schema(path: str) -> list[dict]:
    """Read tables, columns, keys and row counts straight from the database.

    :param path: Path to the database.
    :return: One entry per table.
    """
    conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        names = [
            row[0]
            for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table' AND name NOT LIKE 'sqlite_%' ORDER BY name"
            )
        ]
        tables = []
        for name in names:
            quoted = '"' + name.replace('"', '""') + '"'
            foreign_keys = {
                row[3]: {"table": row[2], "column": row[4]}
                for row in conn.execute(f"PRAGMA foreign_key_list({quoted})")
            }
            columns = [
                {
                    "name": col_name,
                    "type": col_type or "",
                    "notNull": bool(not_null),
                    "primaryKey": pk_index > 0,
                    "references": foreign_keys.get(col_name),
                }
                for _, col_name, col_type, not_null, _, pk_index in conn.execute(f"PRAGMA table_info({quoted})")
            ]
            row_count = conn.execute(f"SELECT COUNT(*) FROM {quoted}").fetchone()[0]
            tables.append({"name": name, "rowCount": row_count, "columns": columns})
        return tables
    finally:
        conn.close()


@blueprint_dataset.route("/api/dataset/info", methods=["GET"])
def dataset_info() -> Response:
    """Describe the downloadable database: file name, size, checksum and schema."""
    try:
        stat = os.stat(DB_PATH)
        return ResponseData(Status.Success, payload={
            "fileName": _download_name(),
            "version": get_version(),
            "sizeBytes": stat.st_size,
            "sha256": _sha256(DB_PATH, stat.st_mtime),
            "tables": _read_schema(DB_PATH),
        }).to_dict()
    except Exception as e:
        return ResponseData(Status.Failure, message=f"failed to read database: {str(e)}").to_dict()


@blueprint_dataset.route("/api/dataset/download", methods=["GET"])
def dataset_download() -> Response:
    """Send the full SQLite database as a file download."""
    # conditional=True gives ETag/Last-Modified and range requests, so resumed
    # or repeated downloads don't resend the whole file
    return send_file(
        DB_PATH,
        mimetype="application/vnd.sqlite3",
        as_attachment=True,
        download_name=_download_name(),
        conditional=True,
    )
