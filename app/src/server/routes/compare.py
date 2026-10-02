# -*- coding: utf-8 -*-

"""Side-by-side comparison of (extended) signatures from the reference database."""

from __future__ import annotations

from functools import lru_cache
from typing import Any

import numpy as np
from flask import Blueprint, jsonify, request
from sqlalchemy import func, or_, select
from sqlalchemy.orm import selectinload

from parasect.core.constants import PROPERTIES
from parasect.database.build_database import (
    AdenylationDomain,
    Protein,
    ProteinDomainAssociation,
    ProteinSynonym,
    Substrate,
    Taxonomy,
)

from .database import get_db
from .network import _clean_taxon

blueprint_compare = Blueprint("compare", __name__)

MAX_ROWS = 300
MAX_FILTER_VALUES = 200

STANDARD_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

TAXONOMY_RANKS = ["domain", "kingdom", "phylum", "cls", "order", "family", "genus", "species", "strain"]


@lru_cache(maxsize=1)
def _similarity_matrix() -> dict[str, Any]:
    """Pairwise residue similarity in [0, 1] from the physicochemical descriptors.

    Each descriptor is z-scored over the 20 standard amino acids, so no single
    one dominates by its units (volume in cubic angstroms would otherwise swamp
    everything). Similarity is then 1 - distance / largest distance between any
    two residues, making identical residues 1 and the most dissimilar pair 0.
    The gap and unknown tokens keep the averaged values the property file gives
    them, so they come out as middling rather than as a match for anything.
    """
    alphabet = "".join(aa for aa in PROPERTIES if len(aa) == 1)
    values = np.array([PROPERTIES[aa] for aa in alphabet], dtype=float)
    standard = np.array([PROPERTIES[aa] for aa in STANDARD_AMINO_ACIDS], dtype=float)
    mean = standard.mean(axis=0)
    std = standard.std(axis=0)
    std[std == 0] = 1.0
    scaled = (values - mean) / std

    distances = np.sqrt(((scaled[:, None, :] - scaled[None, :, :]) ** 2).sum(axis=2))
    similarity = 1.0 - distances / distances.max()
    return {
        "alphabet": alphabet,
        "matrix": np.round(similarity, 3).tolist(),
        "n_properties": int(values.shape[1]),
    }


@blueprint_compare.route("/api/compare/aa_similarity", methods=["GET"])
def aa_similarity():
    """Residue-by-residue similarity matrix used to colour the comparison."""
    return jsonify(_similarity_matrix())


def _values(raw: Any, field: str) -> list[str]:
    """A filter parameter as a de-duplicated list of non-empty strings."""
    if raw is None:
        return []
    if not isinstance(raw, list):
        raw = [raw]
    values = list(dict.fromkeys(str(v).strip() for v in raw if str(v).strip()))
    if len(values) > MAX_FILTER_VALUES:
        raise ValueError(f"at most {MAX_FILTER_VALUES} {field} values at a time")
    return values


def _domain_row(domain: AdenylationDomain) -> dict[str, Any]:
    """Everything the compare page can show for one reference domain."""
    association = domain.proteins[0] if domain.proteins else None
    protein = association.protein if association else None
    taxonomy = protein.taxonomy if protein else None
    return {
        "id": domain.id,
        "name": domain.get_name() or f"domain_{domain.id}",
        "protein": protein.get_name() if protein else "",
        "domain_number": association.domain_number if association else None,
        "signature": domain.signature or "",
        "extended_signature": domain.extended_signature or "",
        "substrates": sorted(s.name for s in domain.substrates),
        "taxonomy": {
            rank: _clean_taxon(getattr(taxonomy, rank, None)) if taxonomy else "unknown"
            for rank in TAXONOMY_RANKS
        },
    }


@blueprint_compare.route("/api/compare/domains", methods=["POST"])
def compare_domains():
    """Reference domains matching any of the given IDs, protein names, substrates or species."""
    body = request.get_json(silent=True) or {}
    try:
        domain_ids = [int(v) for v in (body.get("domain_ids") or [])][:MAX_ROWS]
        proteins = _values(body.get("proteins"), "protein")
        substrates = _values(body.get("substrates"), "substrate")
        species = _values(body.get("species"), "species")
    except (TypeError, ValueError) as e:
        return jsonify({"error": str(e)}), 400

    conditions = []
    if domain_ids:
        conditions.append(AdenylationDomain.id.in_(domain_ids))
    if proteins:
        lowered = [p.lower() for p in proteins]
        conditions.append(AdenylationDomain.proteins.any(ProteinDomainAssociation.protein.has(
            Protein.synonyms.any(func.lower(ProteinSynonym.synonym).in_(lowered)))))
    if substrates:
        lowered = [s.lower() for s in substrates]
        conditions.append(AdenylationDomain.substrates.any(func.lower(Substrate.name).in_(lowered)))
    if species:
        lowered = [s.lower() for s in species]
        conditions.append(AdenylationDomain.proteins.any(ProteinDomainAssociation.protein.has(
            Protein.taxonomy.has(func.lower(Taxonomy.species).in_(lowered)))))
    if not conditions:
        return jsonify({"rows": [], "truncated": False})

    session_generator = get_db()
    session = next(session_generator)
    try:
        query = (
            select(AdenylationDomain)
            .where(or_(*conditions))
            .options(
                selectinload(AdenylationDomain.substrates),
                selectinload(AdenylationDomain.proteins)
                .selectinload(ProteinDomainAssociation.protein)
                .selectinload(Protein.synonyms),
                selectinload(AdenylationDomain.proteins)
                .selectinload(ProteinDomainAssociation.protein)
                .selectinload(Protein.taxonomy),
            )
            .order_by(AdenylationDomain.id)
            .limit(MAX_ROWS + 1)
        )
        domains = list(session.scalars(query).unique())
        truncated = len(domains) > MAX_ROWS
        rows = [_domain_row(d) for d in domains[:MAX_ROWS]]
    finally:
        session_generator.close()

    # explicitly requested IDs first, in the caller's order (a neighbour list
    # arrives ranked, and that ranking is the point)
    position = {domain_id: i for i, domain_id in enumerate(domain_ids)}
    rows.sort(key=lambda r: position.get(r["id"], len(position)))
    return jsonify({"rows": rows, "truncated": truncated})
