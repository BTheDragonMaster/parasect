# -*- coding: utf-8 -*-

"""Sequence-similarity network over the reference PARAS/PARASECT database.

Nodes are adenylation domains from the reference database, related by Hamming
distance between their 34-residue extended signatures. The full domain list
and pairwise distance matrix are computed once per process and cached in
memory: the reference database only changes via redeploys (new entries go
through review as a GitHub PR/issue, see routes/annotation_editor.py), so
there's nothing to invalidate at runtime.

To keep the graph renderable in a browser, the default view groups domains
into clusters (connected components at a chosen Hamming-distance threshold,
the classic "sequence similarity network" technique) rather than showing all
domains at once. Clicking a cluster expands it into its individual members.
"""

from __future__ import annotations

import threading
from typing import Any

import numpy as np
from flask import Blueprint, jsonify, request
from sqlalchemy import select
from sqlalchemy.orm import selectinload

from parasect.database.build_database import AdenylationDomain, Protein, ProteinDomainAssociation

from .database import get_db

blueprint_network = Blueprint("network", __name__)

DEFAULT_THRESHOLD = 5
MAX_THRESHOLD = 20
CLUSTER_EXPAND_KNN = 6  # cap edges-per-node when expanding a large cluster
CLUSTER_EXPAND_FULL_LIMIT = 60  # clusters up to this size get full pairwise edges

_cache_lock = threading.Lock()
_cache: dict[str, Any] | None = None

_cluster_cache_lock = threading.Lock()
_cluster_cache: dict[int, dict[str, Any]] = {}


def _label_counts(label_lists: list[list[str]]) -> dict[str, int]:
    """Count how many domains carry each label (a multi-substrate domain counts once per substrate)."""
    counts: dict[str, int] = {}
    for labels in label_lists:
        for label in set(labels) or {"unknown"}:
            counts[label] = counts.get(label, 0) + 1
    return counts


def _dominant(label_lists: list[list[str]], totals: dict[str, int]) -> str:
    """Return the majority-vote label over a set of domains.

    Every domain casts exactly one vote, split evenly when it carries several
    labels: a handful of reference domains are annotated with more than one
    substrate, and giving each of those a full vote would let one promiscuous
    domain outweigh several specific ones. Ties break on the label's frequency
    across the whole database and then alphabetically, so the same cluster
    always comes out the same colour.
    """
    votes: dict[str, float] = {}
    for labels in label_lists:
        weight = 1.0 / len(labels) if labels else 1.0
        for label in labels or ["unknown"]:
            votes[label] = votes.get(label, 0.0) + weight
    if not votes:
        return "unknown"
    return min(votes.items(), key=lambda kv: (-kv[1], -totals.get(kv[0], 0), kv[0]))[0]


def _clean_taxon(value: str | None) -> str:
    """Normalise the placeholders the database build leaves behind into 'unknown'.

    Some taxonomy rows carry the literal string "None" where lineage lookup
    failed; left alone it shows up in the legend as though it were a real genus.
    """
    if not value or value.strip().lower() in {"none", "na", "n/a", "unclassified"}:
        return "unknown"
    return value


def _build_cache() -> dict[str, Any]:
    """Load every domain's identity/metadata and build the pairwise Hamming distance matrix."""
    session_generator = get_db()
    session = next(session_generator)
    try:
        query = select(AdenylationDomain).options(
            selectinload(AdenylationDomain.synonyms),
            selectinload(AdenylationDomain.substrates),
            selectinload(AdenylationDomain.proteins)
            .selectinload(ProteinDomainAssociation.protein)
            .selectinload(Protein.taxonomy),
        )
        domains = list(session.scalars(query).unique())
    finally:
        session_generator.close()

    ids: list[int] = []
    names: list[str] = []
    substrate_lists: list[list[str]] = []
    genus_list: list[str] = []
    kingdom_list: list[str] = []
    signature_rows: list[list[int]] = []

    for domain in domains:
        sig = domain.extended_signature or ""
        if len(sig) != 34:
            # defensively skip malformed entries rather than let them corrupt the matrix
            continue

        ids.append(domain.id)
        names.append(domain.get_name() or f"domain_{domain.id}")
        substrate_lists.append(sorted({s.name for s in domain.substrates}))

        genus = "unknown"
        kingdom = "unknown"
        if domain.proteins:
            protein = domain.proteins[0].protein
            if protein and protein.taxonomy:
                genus = _clean_taxon(protein.taxonomy.genus)
                kingdom = _clean_taxon(protein.taxonomy.kingdom)
        genus_list.append(genus)
        kingdom_list.append(kingdom)

        signature_rows.append([ord(c) for c in sig])

    signature_matrix = np.array(signature_rows, dtype=np.uint8)
    n = len(ids)

    distances = np.zeros((n, n), dtype=np.uint8)
    chunk = 256
    for start in range(0, n, chunk):
        end = min(start + chunk, n)
        block = (signature_matrix[start:end, None, :] != signature_matrix[None, :, :]).sum(axis=2)
        distances[start:end, :] = block

    return {
        "ids": ids,
        "names": names,
        "substrate_lists": substrate_lists,
        "genus_list": genus_list,
        "kingdom_list": kingdom_list,
        "signature_matrix": signature_matrix,
        "distances": distances,
        "id_to_index": {domain_id: i for i, domain_id in enumerate(ids)},
        # database-wide domain counts per category: the legend uses them to pick
        # its defaults and to break ties in the majority vote
        "totals": {
            "substrate": _label_counts(substrate_lists),
            "genus": _label_counts([[g] for g in genus_list]),
            "kingdom": _label_counts([[k] for k in kingdom_list]),
        },
    }


def _get_cache() -> dict[str, Any]:
    global _cache
    with _cache_lock:
        if _cache is None:
            _cache = _build_cache()
        return _cache


def _connected_components(distances: np.ndarray, threshold: int) -> list[list[int]]:
    """Group indices into connected components where distance <= threshold."""
    n = distances.shape[0]
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    adjacency = distances <= threshold
    for i in range(n):
        neighbours = np.nonzero(adjacency[i, i + 1:])[0] + i + 1
        for j in neighbours:
            union(i, int(j))

    groups: dict[int, list[int]] = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)

    return list(groups.values())


def _get_clusters(threshold: int) -> dict[str, Any]:
    """Compute (and cache) connected-component clusters at a given threshold."""
    with _cluster_cache_lock:
        if threshold in _cluster_cache:
            return _cluster_cache[threshold]

        cache = _get_cache()
        totals = cache["totals"]
        distances = cache["distances"]
        components = _connected_components(distances, threshold)

        clusters = []
        representative_indices = []
        # representative = medoid (member with smallest average distance to the rest)
        for component in components:
            if len(component) == 1:
                representative = component[0]
            else:
                sub = distances[np.ix_(component, component)]
                representative = component[int(np.argmin(sub.sum(axis=1)))]
            representative_indices.append(representative)

            substrates = [cache["substrate_lists"][idx] for idx in component]
            genera = [[cache["genus_list"][idx]] for idx in component]
            kingdoms = [[cache["kingdom_list"][idx]] for idx in component]
            substrate_counts = _label_counts(substrates)

            clusters.append({
                "cluster_id": len(clusters),
                "size": len(component),
                "representative_id": cache["ids"][representative],
                "representative_name": cache["names"][representative],
                "dominant_substrate": _dominant(substrates, totals["substrate"]),
                "dominant_genus": _dominant(genera, totals["genus"]),
                "dominant_kingdom": _dominant(kingdoms, totals["kingdom"]),
                # full composition, so the client can highlight a category that is
                # present in a cluster without winning its majority vote
                "substrate_counts": substrate_counts,
                "genus_counts": _label_counts(genera),
                "kingdom_counts": _label_counts(kingdoms),
                "substrate_diversity": len(substrate_counts),
            })

        # nearest-neighbour edge between each cluster and its single closest other
        # cluster (by representative-to-representative distance), skipped if that
        # distance is far enough that drawing it wouldn't be meaningful
        edges = []
        seen_pairs = set()
        max_edge_distance = threshold + max(5, threshold)
        k = len(clusters)
        if k > 1:
            rep_indices = np.array(representative_indices)
            rep_distances = distances[np.ix_(rep_indices, rep_indices)].astype(np.int16)
            np.fill_diagonal(rep_distances, np.iinfo(np.int16).max)
            for i in range(k):
                j = int(np.argmin(rep_distances[i]))
                dist = int(rep_distances[i, j])
                if dist > max_edge_distance:
                    continue
                pair = (min(i, j), max(i, j))
                if pair in seen_pairs:
                    continue
                seen_pairs.add(pair)
                edges.append({"source": i, "target": j, "distance": dist})

        result = {"clusters": clusters, "edges": edges, "components": components}
        _cluster_cache[threshold] = result
        return result


@blueprint_network.route("/api/network/graph", methods=["GET"])
def get_graph():
    """Return the top-level cluster graph at a given Hamming-distance threshold."""
    try:
        threshold = int(request.args.get("threshold", DEFAULT_THRESHOLD))
    except ValueError:
        return jsonify({"error": "threshold must be an integer"}), 400
    threshold = max(0, min(threshold, MAX_THRESHOLD))

    cache = _get_cache()
    result = _get_clusters(threshold)

    return jsonify({
        "threshold": threshold,
        "total_domains": len(cache["ids"]),
        "clusters": result["clusters"],
        "edges": result["edges"],
    })


@blueprint_network.route("/api/network/cluster/<int:cluster_id>", methods=["GET"])
def expand_cluster(cluster_id: int):
    """Return individual domain nodes/edges for one cluster."""
    try:
        threshold = int(request.args.get("threshold", DEFAULT_THRESHOLD))
    except ValueError:
        return jsonify({"error": "threshold must be an integer"}), 400
    threshold = max(0, min(threshold, MAX_THRESHOLD))

    cache = _get_cache()
    result = _get_clusters(threshold)

    if cluster_id < 0 or cluster_id >= len(result["components"]):
        return jsonify({"error": "cluster not found"}), 404

    member_indices = result["components"][cluster_id]
    distances = cache["distances"]

    nodes = [
        {
            "id": cache["ids"][idx],
            "name": cache["names"][idx],
            "substrates": cache["substrate_lists"][idx],
            # a domain annotated with several substrates has no majority of its own,
            # so the same vote runs over it alone and the tie-break picks the
            # substrate that is most common database-wide
            "dominant_substrate": _dominant([cache["substrate_lists"][idx]], cache["totals"]["substrate"]),
            "genus": cache["genus_list"][idx],
            "kingdom": cache["kingdom_list"][idx],
            "extended_signature": "".join(chr(c) for c in cache["signature_matrix"][idx]),
        }
        for idx in member_indices
    ]

    edges = []
    if len(member_indices) <= CLUSTER_EXPAND_FULL_LIMIT:
        for a in range(len(member_indices)):
            for b in range(a + 1, len(member_indices)):
                dist = int(distances[member_indices[a], member_indices[b]])
                if dist <= threshold:
                    edges.append({
                        "source": cache["ids"][member_indices[a]],
                        "target": cache["ids"][member_indices[b]],
                        "distance": dist,
                    })
    else:
        # large cluster: cap to each node's k nearest neighbours within the cluster
        sub = distances[np.ix_(member_indices, member_indices)]
        seen_pairs = set()
        for a in range(len(member_indices)):
            order = np.argsort(sub[a])
            neighbours = [b for b in order if b != a][:CLUSTER_EXPAND_KNN]
            for b in neighbours:
                pair = (min(a, b), max(a, b))
                if pair in seen_pairs:
                    continue
                seen_pairs.add(pair)
                edges.append({
                    "source": cache["ids"][member_indices[pair[0]]],
                    "target": cache["ids"][member_indices[pair[1]]],
                    "distance": int(sub[pair[0], pair[1]]),
                })

    return jsonify({"cluster_id": cluster_id, "threshold": threshold, "nodes": nodes, "edges": edges})


CATEGORY_FIELDS = ("substrate", "genus", "kingdom")


@blueprint_network.route("/api/network/categories", methods=["GET"])
def get_categories():
    """Every value a colour field takes in the reference database, with domain counts.

    The legend seeds itself with the most frequent categories, but there are ~280
    distinct substrates and ~185 genera in there: picking a rare one needs the
    full list, not just what happens to dominate a cluster in the current view.
    """
    field = (request.args.get("field") or "substrate").strip().lower()
    if field not in CATEGORY_FIELDS:
        return jsonify({"error": f"field must be one of {', '.join(CATEGORY_FIELDS)}"}), 400

    totals = _get_cache()["totals"][field]
    ordered = sorted(totals.items(), key=lambda kv: (-kv[1], kv[0]))
    return jsonify({
        "field": field,
        "categories": [{"label": label, "count": count} for label, count in ordered],
    })


@blueprint_network.route("/api/network/search_names", methods=["GET"])
def search_names():
    """Autocomplete: find domains whose name contains the query substring."""
    q = (request.args.get("q") or "").strip().lower()
    try:
        limit = min(int(request.args.get("limit", 10)), 50)
    except ValueError:
        limit = 10

    if not q:
        return jsonify({"matches": []})

    cache = _get_cache()
    matches = []
    for domain_id, name in zip(cache["ids"], cache["names"]):
        if q in name.lower():
            matches.append({"id": domain_id, "name": name})
            if len(matches) >= limit:
                break

    return jsonify({"matches": matches})


_VALID_AA = set("ACDEFGHIKLMNPQRSTVWYX-")


@blueprint_network.route("/api/network/neighbors", methods=["GET"])
def get_neighbors():
    """Find nearest neighbors of a domain (by ID) or an arbitrary 34-residue signature."""
    try:
        threshold = int(request.args.get("threshold", DEFAULT_THRESHOLD))
    except ValueError:
        return jsonify({"error": "threshold must be an integer"}), 400
    threshold = max(0, min(threshold, MAX_THRESHOLD))

    try:
        k = min(int(request.args.get("k", 10)), 50)
    except ValueError:
        k = 10

    cache = _get_cache()
    domain_id_param = request.args.get("domain_id")
    signature_param = (request.args.get("signature") or "").strip().upper()

    query_index: int | None = None
    query_name = None
    query_signature = None

    if domain_id_param is not None:
        try:
            domain_id = int(domain_id_param)
        except ValueError:
            return jsonify({"error": "domain_id must be an integer"}), 400
        query_index = cache["id_to_index"].get(domain_id)
        if query_index is None:
            return jsonify({"error": "domain not found"}), 404
        query_name = cache["names"][query_index]
        query_signature = "".join(chr(c) for c in cache["signature_matrix"][query_index])
        row_distances = cache["distances"][query_index]
    elif signature_param:
        if len(signature_param) != 34 or not set(signature_param) <= _VALID_AA:
            return jsonify({"error": "signature must be a 34-residue amino acid sequence"}), 400
        query_signature = signature_param
        query_vector = np.array([ord(c) for c in signature_param], dtype=np.uint8)
        row_distances = (cache["signature_matrix"] != query_vector[None, :]).sum(axis=1)
    else:
        return jsonify({"error": "provide either domain_id or signature"}), 400

    order = np.argsort(row_distances)
    neighbours = []
    for idx in order:
        idx = int(idx)
        if idx == query_index:
            continue
        neighbours.append({
            "id": cache["ids"][idx],
            "name": cache["names"][idx],
            "distance": int(row_distances[idx]),
            "substrates": cache["substrate_lists"][idx],
            "genus": cache["genus_list"][idx],
        })
        if len(neighbours) >= k:
            break

    cluster_id = None
    if query_index is not None:
        clusters = _get_clusters(threshold)
        for component_idx, component in enumerate(clusters["components"]):
            if query_index in component:
                cluster_id = component_idx
                break

    return jsonify({
        "query": {"id": int(domain_id_param) if domain_id_param is not None else None,
                   "name": query_name, "signature": query_signature},
        "cluster_id": cluster_id,
        "threshold": threshold,
        "neighbors": neighbours,
    })
