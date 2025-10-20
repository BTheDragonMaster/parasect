# routes/sql.py
# -*- coding: utf-8 -*-
"""
Read-only SQL query endpoints for SQLite.

POST /api/sql         -> paginated JSON for DataGrid
POST /api/sql/export  -> streamed CSV/TSV download
"""

from __future__ import annotations
import os
import re
import sqlite3
import time
from typing import List, Optional

from flask import Blueprint, request, jsonify, abort, Response

try:
    from routes.database import engine as sa_engine  # Optional
except Exception:
    sa_engine = None

blueprint_sql = Blueprint("blueprint_sql", __name__)

DB_PATH = os.getenv("SQLITE_PATH", None)  # can be set via env
MAX_PAGE_SIZE = int(os.getenv("MAX_PAGE_SIZE", "1000"))
MAX_EXPORT_ROWS = int(os.getenv("MAX_EXPORT_ROWS", "100000"))
QUERY_TIMEOUT_SECS = float(os.getenv("QUERY_TIMEOUT_SECS", "30"))

SELECT_LIKE = re.compile(r"^\s*(select|with)\b", re.IGNORECASE)
IDENT = re.compile(r"[A-Za-z_][A-Za-z0-9_]*$")


def _resolve_sqlite_path() -> str:
    """Prefer SQLITE_PATH; else, try to derive from SQLAlchemy engine if it's SQLite; else fallback."""
    if DB_PATH:
        return DB_PATH
    try:
        if sa_engine is not None and sa_engine.url and sa_engine.url.get_backend_name().startswith("sqlite"):
            # this is usually an absolute path; can be ":memory:" too
            return sa_engine.url.database or "records.db"
    except Exception:
        pass
    return "records.db"


def _connect_ro() -> sqlite3.Connection:
    """Open read-only SQLite with a progress handler timeout."""
    db_file = _resolve_sqlite_path()
    conn = sqlite3.connect(f"file:{db_file}?mode=ro", uri=True, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    start = time.time()

    def progress():
        # abort long-running statements (simple protection)
        if time.time() - start > QUERY_TIMEOUT_SECS:
            return 1
        return 0

    conn.set_progress_handler(progress, 10000)
    return conn


def _assert_select_only(sql: str):
    if not SELECT_LIKE.match(sql or ""):
        abort(400, "Only SELECT/WITH queries are allowed.")
    if ";" in sql:
        abort(400, "Semicolons are not allowed.")


def _columns_from_query(conn: sqlite3.Connection, base_query: str) -> List[str]:
    cur = conn.execute(f"SELECT * FROM ({base_query}) AS t LIMIT 0")
    return [d[0] for d in cur.description]


def _total_from_query(conn: sqlite3.Connection, base_query: str) -> int:
    return int(conn.execute(f"SELECT COUNT(*) AS n FROM ({base_query}) AS sub").fetchone()[0])


def _sorted_query(base_query: str, sort_by: Optional[str], sort_dir: Optional[str], valid_cols: List[str]) -> str:
    if not sort_by or not sort_dir:
        return f"SELECT * FROM ({base_query}) AS t"
    # very basic identifier guard + must be a column from the result set
    if not IDENT.fullmatch(sort_by) or sort_by not in valid_cols:
        return f"SELECT * FROM ({base_query}) AS t"
    dir_sql = "ASC" if sort_dir.lower() == "asc" else "DESC"
    return f"SELECT * FROM ({base_query}) AS t ORDER BY \"{sort_by}\" {dir_sql}"


@blueprint_sql.post("/api/sql")
def api_sql():
    """
    Body: { query, page, pageSize, sortBy, sortDir }
    Returns: { columns: [{field, headerName}], rows: [dict], total: int }
    """
    data = request.get_json(force=True, silent=False) or {}
    query = (data.get("query") or "").strip()
    page = int(data.get("page") or 0)
    page_size_req = int(data.get("pageSize") or 25)
    sort_by = data.get("sortBy")
    sort_dir = data.get("sortDir")

    _assert_select_only(query)

    page_size = min(page_size_req, MAX_PAGE_SIZE)
    offset = page * page_size

    try:
        conn = _connect_ro()
        cols = _columns_from_query(conn, query)
        total = _total_from_query(conn, query)
        sql_sorted = _sorted_query(query, sort_by, sort_dir, cols)
        rows = conn.execute(f"{sql_sorted} LIMIT ? OFFSET ?", (page_size, offset)).fetchall()
    except sqlite3.Error as e:
        abort(400, f"Query failed: {e}")
    finally:
        try:
            conn.close()
        except Exception:
            pass

    return jsonify({
        "columns": [{"field": c, "headerName": c} for c in cols],
        "rows": [{k: r[k] for k in r.keys()} for r in rows],
        "total": int(total),
    })


@blueprint_sql.post("/api/sql/export")
def api_sql_export():
    """
    Server-side export.
    Body: { query, sortBy, sortDir, format: 'csv'|'tsv' }
    Returns: streamed CSV/TSV with up to MAX_EXPORT_ROWS rows.
    """
    data = request.get_json(force=True, silent=False) or {}
    query = (data.get("query") or "").strip()
    sort_by = data.get("sortBy")
    sort_dir = data.get("sortDir")
    fmt = (data.get("format") or "csv").lower()
    delim = "," if fmt == "csv" else "\t"
    mime = "text/csv" if fmt == "csv" else "text/tab-separated-values"

    _assert_select_only(query)

    try:
        conn = _connect_ro()
        cols = _columns_from_query(conn, query)
        sql_sorted = _sorted_query(query, sort_by, sort_dir, cols)
        cur = conn.execute(f"{sql_sorted} LIMIT ?", (MAX_EXPORT_ROWS,))
        colnames = cols

        def generate():
            # header
            yield delim.join(colnames) + "\n"
            for row in cur:
                vals = []
                for c in colnames:
                    v = row[c]
                    if v is None:
                        vals.append("")
                    else:
                        s = str(v).replace("\r", " ").replace("\n", " ")
                        if delim == "," and ("," in s or '"' in s):
                            s = '"' + s.replace('"', '""') + '"'
                        vals.append(s)
                yield delim.join(vals) + "\n"

        return Response(
            generate(),
            mimetype=mime,
            headers={"Content-Disposition": f'attachment; filename="export.{fmt}"'},
        )
    except sqlite3.Error as e:
        abort(400, f"Export failed: {e}")
    finally:
        try:
            conn.close()
        except Exception:
            pass
