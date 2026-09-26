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


def _columns_from_query(conn: sqlite3.Connection, base_query: str, params=()) -> list[str]:
    cur = conn.execute(f"SELECT * FROM ({base_query}) AS t LIMIT 0", params)
    return [d[0] for d in cur.description]


def _total_from_query(conn: sqlite3.Connection, base_query: str, params=()) -> int:
    return int(conn.execute(f"SELECT COUNT(*) AS n FROM ({base_query}) AS sub", params).fetchone()[0])


def _sorted_query(base_query: str, sort_by: str | None, sort_dir: str | None, valid_cols: list[str]) -> str:
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


########################################################################################################################
#
# Safe preset queries: the query text is always a fixed, hard-coded SQL template chosen from PRESETS below and
# never built from client-supplied text. User input only ever fills bound SQL parameters (:name placeholders),
# which SQLite binds as data, never as SQL syntax. So unlike /api/sql (which runs whatever SELECT text a
# caller sends), this endpoint can't be used to inject SQL regardless of what a user types into a filter field.
#
########################################################################################################################


MAX_FILTER_VALUES = 200


def _as_values(raw, field: str) -> list[str]:
    """Normalise a filter parameter into a list of non-empty strings.

    Filters accept either a single value or a list, so an older client sending
    a bare string keeps working alongside the multi-select UI.
    """
    if raw is None:
        values = []
    elif isinstance(raw, (list, tuple)):
        values = [str(v).strip() for v in raw]
    else:
        values = [str(raw).strip()]
    values = [v for v in values if v]
    if not values:
        raise ValueError(f"{field} is required")
    if len(values) > MAX_FILTER_VALUES:
        raise ValueError(f"at most {MAX_FILTER_VALUES} {field} values at a time")
    # de-duplicate, keeping the caller's order
    return list(dict.fromkeys(values))


def _in_clause(values: list[str], prefix: str):
    """Build an IN (...) list of bound placeholders.

    The placeholder NAMES are generated here (:sub0, :sub1, ...) and only the
    values are bound, so a multi-value filter is exactly as injection-proof as
    the single-value one. No user text ever reaches the SQL text.
    """
    keys = [f"{prefix}{i}" for i in range(len(values))]
    return ", ".join(f":{k}" for k in keys), dict(zip(keys, values))


def _build_by_substrate(params: dict):
    names = _as_values(params.get("substrate_name"), "substrate_name")
    placeholders, bound = _in_clause(names, "sub")
    sql = """
        SELECT
          ad.id                          AS domain_id,
          p.id                           AS protein_id,
          group_concat(ps.synonym, ', ') AS protein_synonyms,
          pda.domain_number,
          ad.signature,
          ad.extended_signature,
          s.name                         AS substrate_name,
          s.smiles                       AS substrate_smiles
        FROM substrate_domain_association sda
        JOIN substrate                    s   ON s.name = sda.substrate_name
        JOIN adenylation_domain           ad  ON ad.id = sda.domain_id
        JOIN protein_domain_association   pda ON pda.domain_id = ad.id
        JOIN protein                      p   ON p.id = pda.protein_id
        LEFT JOIN protein_synonym         ps  ON ps.protein_id = p.id
        WHERE s.name COLLATE NOCASE IN (PLACEHOLDERS)
        GROUP BY ad.id, p.id, pda.domain_number, ad.signature, ad.extended_signature, s.name, s.smiles
        ORDER BY s.name, p.id, pda.domain_number
    """.replace("PLACEHOLDERS", placeholders)
    return sql, bound


def _build_by_protein_id(params: dict):
    protein_ids = _as_values(params.get("protein_id"), "protein_id")
    placeholders, bound = _in_clause(protein_ids, "prot")
    sql = """
        SELECT
          ad.id                     AS domain_id,
          p.id                      AS protein_id,
          ps.synonym                AS protein_synonym,
          pda.domain_number,
          ad.signature,
          ad.extended_signature,
          s.name                    AS substrate_name,
          s.smiles                  AS substrate_smiles
        FROM protein_synonym              ps
        JOIN protein                      p    ON p.id = ps.protein_id
        JOIN protein_domain_association   pda  ON pda.protein_id = p.id
        JOIN adenylation_domain           ad   ON ad.id = pda.domain_id
        LEFT JOIN substrate_domain_association sda ON sda.domain_id = ad.id
        LEFT JOIN substrate               s    ON s.name = sda.substrate_name
        WHERE ps.synonym COLLATE NOCASE IN (PLACEHOLDERS)
        ORDER BY ps.synonym, pda.domain_number, s.name
    """.replace("PLACEHOLDERS", placeholders)
    return sql, bound


def _build_by_species(params: dict):
    species_names = _as_values(params.get("species"), "species")
    placeholders, bound = _in_clause(species_names, "sp")
    sql = """
        SELECT
          ad.id                     AS domain_id,
          p.id                      AS protein_id,
          t.species                 AS species,
          pda.domain_number,
          ad.signature,
          ad.extended_signature,
          s.name                    AS substrate_name,
          s.smiles                  AS substrate_smiles
        FROM taxonomy                     t
        JOIN protein                      p    ON p.taxonomy_id = t.id
        JOIN protein_domain_association   pda  ON pda.protein_id = p.id
        JOIN adenylation_domain           ad   ON ad.id = pda.domain_id
        LEFT JOIN substrate_domain_association sda ON sda.domain_id = ad.id
        LEFT JOIN substrate               s    ON s.name = sda.substrate_name
        WHERE t.species COLLATE NOCASE IN (PLACEHOLDERS)
        ORDER BY t.species, pda.domain_number, s.name
    """.replace("PLACEHOLDERS", placeholders)
    return sql, bound


def _build_by_signature(params: dict):
    signature = (params.get("signature") or "").strip().upper()
    if not signature:
        raise ValueError("signature is required")
    try:
        max_distance = int(params.get("max_distance", 3))
    except (TypeError, ValueError):
        raise ValueError("max_distance must be an integer")
    max_distance = max(0, min(10, max_distance))
    padded = signature[:10].ljust(10, "-")  # '-' is the gap/wildcard token

    # 10 position-by-position comparisons against the query signature; :padded is
    # bound once and safely reused across all 10 substr() calls
    comparisons = " + ".join(
        f"(substr(:padded,{i + 1},1) <> substr(UPPER(substr(ad.signature || '----------',1,10)),{i + 1},1))"
        for i in range(10)
    )
    sql = f"""
        WITH d AS (
            SELECT ad.id, ad.signature, ad.extended_signature, ({comparisons}) AS hamming
            FROM adenylation_domain ad
        )
        SELECT
          d.id                       AS domain_id,
          p.id                       AS protein_id,
          group_concat(ps.synonym, ', ') AS protein_synonyms,
          pda.domain_number,
          ad.signature,
          ad.extended_signature,
          s.name                     AS substrate_name,
          s.smiles                   AS substrate_smiles,
          d.hamming
        FROM d
        JOIN adenylation_domain           ad  ON ad.id = d.id
        JOIN protein_domain_association   pda ON pda.domain_id = ad.id
        JOIN protein                      p   ON p.id = pda.protein_id
        LEFT JOIN substrate_domain_association sda ON sda.domain_id = ad.id
        LEFT JOIN substrate               s   ON s.name = sda.substrate_name
        LEFT JOIN protein_synonym         ps  ON ps.protein_id = p.id
        WHERE d.hamming <= :max_distance
        GROUP BY d.id, p.id, pda.domain_number, ad.signature, ad.extended_signature, s.name, s.smiles, d.hamming
        ORDER BY pda.domain_number, s.name
    """
    return sql, {"padded": padded, "max_distance": max_distance}


PRESETS = {
    "substrate": _build_by_substrate,
    "proteinId": _build_by_protein_id,
    "species": _build_by_species,
    "signature": _build_by_signature,
}


# Option lists for the filter dropdowns. Like the preset queries these are fixed
# SQL text: the only user input is a search string, bound as a LIKE pattern.
_OPTION_QUERIES = {
    "substrate": """
        SELECT s.name AS value, COUNT(DISTINCT sda.domain_id) AS count
        FROM substrate s
        LEFT JOIN substrate_domain_association sda ON sda.substrate_name = s.name
        WHERE :q = '' OR s.name LIKE :like ESCAPE '\\'
        GROUP BY s.name
        ORDER BY count DESC, s.name
        LIMIT :limit
    """,
    "protein": """
        SELECT ps.synonym AS value, COUNT(DISTINCT pda.domain_id) AS count
        FROM protein_synonym ps
        JOIN protein_domain_association pda ON pda.protein_id = ps.protein_id
        WHERE :q = '' OR ps.synonym LIKE :like ESCAPE '\\'
        GROUP BY ps.synonym
        ORDER BY count DESC, ps.synonym
        LIMIT :limit
    """,
    "species": """
        SELECT t.species AS value, COUNT(DISTINCT pda.domain_id) AS count
        FROM taxonomy t
        JOIN protein p ON p.taxonomy_id = t.id
        JOIN protein_domain_association pda ON pda.protein_id = p.id
        WHERE :q = '' OR t.species LIKE :like ESCAPE '\\'
        GROUP BY t.species
        ORDER BY count DESC, t.species
        LIMIT :limit
    """,
}


@blueprint_sql.get("/api/sql/options")
def api_sql_options():
    """Searchable option list for a filter dropdown.

    Body-less GET: ?field=substrate|protein|species&q=<search>&limit=<n>.

    The dropdowns load through here rather than pulling every value up front:
    there are already ~1,900 distinct protein identifiers, and that list only
    grows as the reference database does.
    """
    field = (request.args.get("field") or "").strip().lower()
    sql = _OPTION_QUERIES.get(field)
    if sql is None:
        abort(400, f"Unknown field '{field}'. Valid fields: {sorted(_OPTION_QUERIES)}")

    q = (request.args.get("q") or "").strip()
    try:
        limit = min(max(int(request.args.get("limit", 50)), 1), 200)
    except ValueError:
        limit = 50

    # neutralise LIKE wildcards so a search for "P_1" matches that literally
    escaped = q.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")

    try:
        conn = _connect_ro()
        rows = conn.execute(sql, {"q": q, "like": f"%{escaped}%", "limit": limit}).fetchall()
    except sqlite3.Error as e:
        abort(400, f"Option lookup failed: {e}")
    finally:
        try:
            conn.close()
        except Exception:
            pass

    return jsonify({
        "field": field,
        "options": [{"value": r["value"], "count": int(r["count"])} for r in rows if r["value"]],
    })


# Aggregate counts behind the summary chart on the query page. Same deal as the
# option lists above: the SQL text is fixed and there is no user input at all,
# the client only picks one of these views by key.

_ANNOTATED_DOMAINS = """
    FROM adenylation_domain           ad
    JOIN substrate_domain_association sda ON sda.domain_id = ad.id          -- annotated only
    JOIN protein_domain_association   pda ON pda.domain_id = ad.id
    JOIN protein                      p   ON p.id = pda.protein_id
    JOIN taxonomy                     t   ON t.id = p.taxonomy_id
"""


def _by_taxonomic_rank(column: str) -> str:
    """Build the per-rank count query for one taxonomy column.

    `column` never comes from a request as the only call sites are the literal
    names written into _STATS_QUERIES below, so nothing user-supplied can reach
    the SQL text here.

    Ranks that were never assigned arrive as the literal string 'None' (that is
    what the loader writes for a missing rank, not SQL NULL), which reads as a
    real taxon name on an axis. They are folded into one 'unclassified' bucket
    so the chart doesn't claim a genus called None.

    :param column: name of the column on `taxonomy` to group by.
    :type column: str
    :returns: SQL selecting label/value rows, largest first.
    :rtype: str
    """
    return f"""
        SELECT
          CASE
            WHEN t.{column} IS NULL OR trim(t.{column}) = ''
              OR lower(t.{column}) IN ('none', 'unknown', 'na', 'n/a')
            THEN 'unclassified'
            ELSE t.{column}
          END                   AS label,
          COUNT(DISTINCT ad.id) AS value
        {_ANNOTATED_DOMAINS}
        GROUP BY label
        ORDER BY value DESC, label
    """


_STATS_QUERIES = {
    "domain": _by_taxonomic_rank("domain"),
    "phylum": _by_taxonomic_rank("phylum"),
    "genus": _by_taxonomic_rank("genus"),
    "species": _by_taxonomic_rank("species"),
    "substrate": """
        SELECT
          s.name                        AS label,
          COUNT(DISTINCT sda.domain_id) AS value
        FROM substrate_domain_association sda
        JOIN substrate                    s ON s.name = sda.substrate_name
        GROUP BY s.name
        ORDER BY value DESC, label
    """,
}


@blueprint_sql.get("/api/sql/stats")
def api_sql_stats():
    """Category counts for the summary chart on the query page.

    Body-less GET: ?view=domain|phylum|genus|species|substrate.

    The whole grouping comes back in one response rather than a server-side
    top-N: the widest of these views is ~520 rows, small enough that letting the
    client re-slice how many bars it shows beats another round trip per view.
    """
    view = (request.args.get("view") or "").strip().lower()
    sql = _STATS_QUERIES.get(view)
    if sql is None:
        abort(400, f"Unknown view '{view}'. Valid views: {sorted(_STATS_QUERIES)}")

    try:
        conn = _connect_ro()
        rows = conn.execute(sql).fetchall()
    except sqlite3.Error as e:
        abort(400, f"Stats lookup failed: {e}")
    finally:
        try:
            conn.close()
        except Exception:
            pass

    return jsonify({
        "view": view,
        "rows": [{"label": r["label"], "value": int(r["value"])} for r in rows],
    })


@blueprint_sql.post("/api/sql/preset")
def api_sql_preset():
    """
    Body: { preset: str, params: {...}, page, pageSize, sortBy, sortDir }
    Returns: same shape as /api/sql -> { columns, rows, total }
    """
    data = request.get_json(force=True, silent=False) or {}
    preset_key = data.get("preset")
    params = data.get("params") or {}
    page = int(data.get("page") or 0)
    page_size = min(int(data.get("pageSize") or 25), MAX_PAGE_SIZE)
    sort_by = data.get("sortBy")
    sort_dir = data.get("sortDir")

    builder = PRESETS.get(preset_key)
    if not builder:
        abort(400, f"Unknown preset '{preset_key}'. Valid presets: {sorted(PRESETS)}")

    try:
        base_query, bound_params = builder(params)
    except ValueError as e:
        abort(400, str(e))

    offset = page * page_size

    try:
        conn = _connect_ro()
        cols = _columns_from_query(conn, base_query, bound_params)
        total = _total_from_query(conn, base_query, bound_params)
        sql_sorted = _sorted_query(base_query, sort_by, sort_dir, cols)
        rows = conn.execute(
            f"{sql_sorted} LIMIT :__limit OFFSET :__offset",
            {**bound_params, "__limit": page_size, "__offset": offset},
        ).fetchall()
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
