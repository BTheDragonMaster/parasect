#!/usr/bin/env python3
"""
fetch_ncbi_taxonomy.py

Query NCBI Taxonomy for a list of names and produce a tab-separated TSV with:
supplied_name, matched_name, tax_id, rank, kingdom, phylum, class, order, family, genus, species, notes

Requirements:
 - biopython
 - pandas
 - tqdm

Install:
 pip install biopython pandas tqdm

Usage:
 export NCBI_EMAIL="your@email"
 export NCBI_API_KEY="your_ncbi_api_key"   # optional but recommended
 python fetch_ncbi_taxonomy.py --input names.txt --output taxonomy.tsv --batch 4

Notes:
 - NCBI limits requests per second; this script sleeps automatically to respect rate limits.
 - Ambiguous matches are flagged in the notes column.
"""
import os
import sys
import time
import argparse
from Bio import Entrez
import pandas as pd
from tqdm import tqdm

# Configure Entrez
ENTREZ_EMAIL = "barbara.r.terlouw@gmail.com"

Entrez.email = ENTREZ_EMAIL

# Safe delay (NCBI suggests up to ~3 requests/sec with an API key; we'll be conservative)
DELAY = 0.6


def esearch_name(name):
    """Search taxonomy names; return list of taxids (strings) in order of relevance"""
    # Use scientific name search and fallback to general text search
    query = f'"{name}"[Scientific Name] OR "{name}"[All Names]'
    handle = Entrez.esearch(db="taxonomy", term=query, retmode="xml")
    result = Entrez.read(handle)
    handle.close()
    ids = result.get("IdList", [])
    return ids


def efetch_taxon(taxid):
    """Fetch taxonomy record for a single taxid"""
    handle = Entrez.efetch(db="taxonomy", id=str(taxid), retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    if not records:
        return None
    rec = records[0]
    return rec

def parse_lineage(rec):
    """Return dict with kingdom..species from Entrez taxonomy record"""
    out = dict.fromkeys(["kingdom","phylum","class","order","family","genus","species"], "")
    # rec has 'LineageEx' list of dicts with 'Rank' and 'ScientificName'
    lineage = rec.get("LineageEx", [])
    for node in lineage:
        r = node.get("Rank","").lower()
        n = node.get("ScientificName","")
        if r == "superkingdom" or r == "kingdom":
            out["kingdom"] = n
        elif r == "phylum":
            out["phylum"] = n
        elif r == "class":
            out["class"] = n
        elif r == "order":
            out["order"] = n
        elif r == "family":
            out["family"] = n
        elif r == "genus":
            out["genus"] = n
    # species: use the record's ScientificName if rank==species, otherwise try to parse
    if rec.get("Rank","").lower() == "species":
        out["species"] = rec.get("ScientificName","")
    else:
        # attempt to find a species node in LineageEx
        for node in reversed(lineage):
            if node.get("Rank","").lower() == "species":
                out["species"] = node.get("ScientificName","")
                break
    return out

def best_match_for_name(name):
    """
    For a supplied name, attempt to find the best NCBI Taxonomy match.
    Strategy:
     - esearch by scientific name / all names
     - if single hit -> accept
     - if multiple hits -> fetch top N (default 5) and pick:
         * exact case-insensitive scientific name match first
         * exact match among synonyms / other names
         * otherwise choose top hit but flag as 'ambiguous'
    """
    notes = ""
    ids = esearch_name(name)
    time.sleep(DELAY)
    if not ids:
        # No hits: try a more permissive search (without quotes)
        handle = Entrez.esearch(db="taxonomy", term=name, retmode="xml")
        res = Entrez.read(handle)
        handle.close()
        ids = res.get("IdList", [])
        time.sleep(DELAY)
        notes += "no initial hits; used permissive search. "
    if not ids:
        return {"supplied": name, "matched": "", "taxid": "", "rank":"", "lineage":{}, "notes": "no match found"}
    # Fetch top candidates (up to 5)
    top_ids = ids[:5]
    records = []
    for tid in top_ids:
        rec = efetch_taxon(tid)
        time.sleep(DELAY)
        if rec:
            records.append(rec)
    # Attempt exact scientific name match (case-insensitive)
    lower_name = name.strip().lower()
    for rec in records:
        if rec.get("ScientificName","").lower() == lower_name:
            lineage = parse_lineage(rec)
            return {"supplied": name, "matched": rec.get("ScientificName",""), "taxid": rec.get("TaxId",""),
                    "rank": rec.get("Rank",""), "lineage": lineage, "notes": notes + "exact scientific name match"}
    # Check synonyms and other name lists (OtherNames may include synonyms)
    for rec in records:
        other = []
        on = rec.get("OtherNames", {})
        for key in ("Synonym","EquivalentName","GenbankCommonName","BlastName","GenbankCommonName"):
            vals = on.get(key, [])
            if isinstance(vals, str):
                other.append(vals)
            elif isinstance(vals, list):
                other.extend(vals)
        other = [x.lower() for x in other if x]
        if lower_name in other:
            lineage = parse_lineage(rec)
            return {"supplied": name, "matched": rec.get("ScientificName",""), "taxid": rec.get("TaxId",""),
                    "rank": rec.get("Rank",""), "lineage": lineage, "notes": notes + "matched via synonym/other name"}
    # Otherwise choose top record but mark ambiguous
    rec = records[0]
    lineage = parse_lineage(rec)
    notes += f"ambiguous: {len(records)} candidate(s) - top hit used"
    return {"supplied": name, "matched": rec.get("ScientificName",""), "taxid": rec.get("TaxId",""),
            "rank": rec.get("Rank",""), "lineage": lineage, "notes": notes}

def process_list(names, sleep=DELAY):
    rows = []
    for name in tqdm(names, desc="Processing names"):
        name = name.strip()
        if not name:
            continue
        try:
            res = best_match_for_name(name)
        except Exception as e:
            res = {"supplied": name, "matched": "", "taxid": "", "rank":"", "lineage":{}, "notes": f"error: {e}"}
        line = {
            "supplied_name": res["supplied"],
            "matched_name": res.get("matched",""),
            "tax_id": res.get("taxid",""),
            "rank": res.get("rank",""),
            "kingdom": res.get("lineage",{}).get("kingdom",""),
            "phylum": res.get("lineage",{}).get("phylum",""),
            "class": res.get("lineage",{}).get("class",""),
            "order": res.get("lineage",{}).get("order",""),
            "family": res.get("lineage",{}).get("family",""),
            "genus": res.get("lineage",{}).get("genus",""),
            "species": res.get("lineage",{}).get("species",""),
            "notes": res.get("notes","")
        }
        rows.append(line)
    return pd.DataFrame(rows)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", "-o", default="taxonomy.tsv", help="output TSV file")
    parser.add_argument("--batch", type=int, default=4, help="number of parallel queries (not used currently)")
    args = parser.parse_args()

    names = list({'Bacillaceae', 'Eubacteriales', 'Archangiaceae', 'Nostoc', 'Actinomycetales', 'Hyphomicrobiales', 'Euteleostei', 'Marasmiineae', 'Klebsiella', 'Agaricales', 'Brachycera', 'Diptera', 'Staphylococcus', 'Pyricularia', 'Janthinobacterium', 'Rhizobium/Agrobacterium group', 'Ustilaginomycotina', 'Pyxidicoccus', 'Hypocreales', 'Chloroflexota', 'Aeromonadales', 'Lacticaseibacillus', 'Thermothelomyces', 'Pseudomonadati', 'Didymosphaeriaceae', 'Moorena', 'Myoxocephalus', 'Pseudomonadales', 'Xenorhabdus', 'Jahnella', 'Neonothopanus', 'Herpetosiphonales', 'Trichoderma', 'Cuspidothrix', 'Pseudovibrio', 'Pezizomycotina', 'Polyangiaceae', 'Eurotiomyceta', 'Polyangia', 'Chaetomiaceae', 'Alphaproteobacteria', 'Sorangium', 'Marinactinospora', 'Nodularia', 'Ophiocordycipitaceae', 'Myxococcota', 'Chromobacterium', 'Magnaporthales', 'Yersinia', 'Sordariomycetes', 'Mortierellales', 'Chaetothyriomycetidae', 'Xylariomycetidae', 'Aneurinibacillus group', 'Herpetosiphon', 'Ruminococcaceae', 'Eukaryota', 'Pyrenochaetopsis', 'Aspergillaceae', 'Nonomuraea', 'Cladobotryum', 'Chromobacteriaceae', 'Ephydroidea', 'Coniochaetales', 'Pseudonocardiales', 'Delftia', 'Tistrella', 'Aspergillus subgen. Nidulantes', 'Talaromyces sect. Islandici', 'Neisseriales', 'Oxalobacteraceae', 'Suillineae', 'Pyrenophora', 'Pseudomyxococcus', 'Rhizobiales', 'Listeriaceae', 'OSLEUM clade', 'Goodfellowiella', 'Suillaceae', 'Claviceps', 'Rothia', 'Saccharospirillaceae', 'Neopterygii', 'Nannocystales', 'Eurotiales', 'Aeromonas', 'Brucella', 'Penicillium', 'Paraburkholderia', 'Eurotiomycetes', 'Euteleostomi', 'Schizosaccharomycetales', 'Mucoromycota', 'Melittangium', 'Serpulaceae', 'Clostridia', 'Hypocreales incertae sedis', 'Achaetomiella', 'Fusarium vanettenii', 'Phoma', 'Lysobacter', 'Malbranchea', 'Brucellaceae', 'Moraxellaceae', 'Thalassospira', 'Aspergillus subgen. Fumigati', 'Archangium', 'Okeania', 'Morganellaceae', 'Oceanospirillales', 'Dipodascomycetes', 'Paenibacillaceae', 'C4D1M', 'Agaricineae', 'Drosophilidae', 'Holometabola', 'Scedosporium', 'Parastagonospora', 'sordariomyceta', 'Microcoleaceae', 'Xanthomonadaceae', 'Botryosphaeriales', 'Fusarium solani species complex', 'Hapsidospora', 'Tapinellineae', 'Anabaena', 'Kaarinaea lacus', 'Kaarinaea', 'Pseudoalteromonadaceae', 'Micromonospora', 'Cladosporiaceae', 'Taphrinomycotina', 'Heterobasidion', 'Streptomyces aurantiacus group', 'Oscillospiraceae', 'Aquimarina', 'Cyanobacteriota/Melainabacteria group', 'Paraburkholderia graminis', 'Aphanizomenonaceae', 'Coniochaeta', 'Neoteleostei', 'Rhizobium', 'Thalassospiraceae', 'Helotiales', 'Chloroflexaceae', 'Photorhabdus', 'Periglandula', 'Candidatus Profftella', 'Pleosporaceae', 'Chaetothyriaceae', 'Insecta', 'Cladosporium', 'Nodulariaceae', 'Bacillales', 'Glarea', 'Pleosporomycetidae', 'Polyporales', 'Staphylococcaceae', 'Kibdelosporangium', 'Fungi incertae sedis', 'Eurotiomycetidae', 'Physalacriaceae', 'Candidatus Synoicihabitans', 'Bacillati', 'Polyangiales', 'Wickerhamiella', 'Russulales', 'Orbiliales', 'Glomerellales', 'Pseudomonadaceae', 'Candidatus Tectimicrobiota', 'Discosia', 'Pyrenochaetopsidaceae', 'Malbrancheaceae', 'Lactococcus', 'Amycolatopsis', 'Chloroflexi', 'Teredinibacter', 'Ustilaginaceae', 'Chroococcales', 'Enterobacteriaceae', 'Microcystis', 'Phyllosticta', 'Chaetothyriales', 'Thermomonosporaceae', 'Suillus', 'Vertebrata', 'Fusarium oxysporum species complex', 'Colletotrichum destructivum species complex', 'Burkholderiaceae', 'Cystobacterineae', 'Neurospora', 'Oscillatoriales', 'Agaricomycotina', 'Candidatus Entotheonellia', 'Kitasatosporales', 'Listonella', 'Kallichroma', 'Flavobacteriaceae', 'Pseudoalteromonas', 'Stachybotrys', 'Stappiaceae', 'Percomorpha', 'Schizosaccharomyces', 'Pseudopithomyces', 'Trichocomaceae', 'Sordariomycetidae', 'Mariannaea', 'Ruminococcus', 'Schizosaccharomycetaceae', 'Lysobacterales', 'Colletotrichum', 'Psathyrellaceae', 'Oceanobacillus', 'Enterobacterales', 'Schizosaccharomycetes', 'Stachybotryaceae', 'Mycobacterium', 'Drosophila', 'Orbilia', 'Stigmatella', 'Lecanicillium', 'Sporocadaceae', 'Nostocaceae', 'Acremonium', 'Nostocales', 'Streptomyces albidoflavus group', 'Paenibacillus', 'Actinomycetes', 'Actinopterygii', 'Telluria group', 'Leptosphaeriaceae', 'unclassified sequences', 'Rhodococcus', 'Omphalotaceae', 'Teleostei', 'Fusarium fujikuroi species complex', 'Muscomorpha', 'Xenoacremonium', 'Verrucomicrobiota', 'Orbiliomycetes', 'Hyalangium', 'Phaeosphaeriaceae', 'Planktothrix', 'Cordycipitaceae', 'dothideomyceta', 'Escovopsis', 'Oscillatoriaceae', 'Onygenales', 'Rhizobiaceae', 'Curvularia', 'Penicillium chrysogenum species complex', 'Mortierellomycotina', 'Nostoc cyanobionts', 'Chondromyces', 'Gammaproteobacteria', 'artificial sequences', 'Saccharomycetes', 'Corynebacterineae', 'Moraxellales', 'Micromonosporales', 'Actinobacteridae', 'Gliocladium', 'Streptoalloteichus', 'Micromonosporaceae', 'Mycobacteriales', 'Bacteroidota', 'Trametes', 'Flavobacteriales', 'Pleosporineae', 'Chaetomium', 'Arthropoda', 'Alteromonadales', 'Acanthomorpha', 'Botryosphaeriaceae', 'Mortierellaceae', 'Mycobacteroides', 'Coniophorineae', 'Clostridiales', 'Ustilaginomycetes', 'Kitasatospora', 'Cottidae', 'Nannocystaceae', 'Betaproteobacteria', 'Hexapoda', 'Agaricomycetes', 'Halomonas', 'Rhodospirillales', 'Nitrospinota/Tectimicrobiota group', 'Chordata', 'Tolypocladium', 'Dothideomycetes incertae sedis', 'Bacillus amyloliquefaciens group', 'Cottoidei', 'Chloroflexia', 'Saccharomycetaceae', 'Streptosporangiales', 'Saccharomonospora', 'Candidatus Entotheonella', 'Deltaproteobacteria', 'Fungi', 'Boletales', 'Coleofasciculales', 'Scorpaeniformes', 'Polyangium', 'Heterobasidion annosum species complex', 'Hypocreomycetidae', 'Massilia', 'Agaricomycetidae', 'Streptomyces rochei group', 'Ustilaginales', 'Nocardiaceae', 'Fusarium', 'Amphichorda', 'Myxococcales', 'Acinetobacter calcoaceticus/baumannii complex', 'Aspergillus', 'Yersiniaceae', 'Coprinopsis', 'Streptomycetaceae', 'Ralstonia', 'Lecanoromycetes', 'Sordariaceae', 'Aspergillus subgen. Circumdati', 'Leotiomycetes', 'Xanthomonas', 'Mycobacteriaceae', 'Agrobacterium', 'Fusarium incarnatum-equiseti species complex', 'Phyllostictaceae', 'Epichloe', 'Hapalosiphonaceae', 'Glomerellaceae', 'Serpula', 'Fusarium heterosporum species complex', 'Bacillus cereus group', 'Sandaracinaceae', 'Craniata', 'Pseudomonas', 'Helotiaceae', 'Cellvibrionales', 'Saccharomyces', 'Mycosarcoma', 'Microascales', 'Trichomonascaceae', 'Sordariales', 'Agrobacterium tumefaciens complex', 'Caldimonas', 'Alternaria alternata complex', 'Pyriculariaceae', 'Dothideomycetes', 'Acanthopterygii', 'Streptomyces violaceusniger group', 'Bondarzewiaceae', 'Metazoa', 'Streptomyces', 'Coleofasciculaceae', 'Pterygota', 'Saccharothrix', 'Kutzneria', 'Talaromyces', 'Cupriavidus', 'Micrococcaceae', 'Dothideomycetidae', 'Pestalotiopsis', 'Fischerella', 'Streptosporangium', 'Actinobacteria', 'Salmonella', 'vectors', 'Serratia', 'Ralstonia solanacearum species complex', 'Tapinellaceae', 'Mortierella', 'Cladosporiales', 'Fusarium tricinctum species complex', 'Firmicutes', 'Vibrionales', 'Cyanophyceae', 'Clavicipitaceae', 'Aeromonadaceae', 'Opisthokonta', 'Burkholderia', 'Ostropomycetidae', 'Rhodococcus erythropolis group', 'Gynuella', 'Mycobacterium tuberculosis complex', 'Geminicoccaceae', 'Sorangiineae', 'Actinosynnema', 'Ostropales', 'Proteobacteria', 'Escherichia', 'Lactiplantibacillus', 'Cylindrospermopsis', 'Cyanobacteriota', 'Nannocystis', 'Opitutales', 'Cordyceps', 'Opitutaceae', 'Lyngbya', 'Salinispora', 'Mortierellomycetes', 'Amphisphaeriales', 'Phormidium', 'Vibrio', 'Dikarya', 'Alternaria', 'saccharomyceta', 'Burkholderia cepacia complex', 'Bacillus', 'Streptomycetales', 'Herpetosiphonaceae', 'Sophophora', 'Niallia', 'Azotobacter', 'Streptococcaceae', 'Stictidaceae', 'Mycolicibacterium', 'Cylindrospermum', 'Trichormus', 'Saccotheciaceae', 'Cystobacter', 'leotiomyceta', 'Coniochaetaceae', 'Corynebacteriales', 'Flavobacteriia', 'Bacteria', 'Pseudomonas syringae', 'Aneurinibacillus', 'Vibrionaceae', 'Micrococcales', 'Ruminiclostridium', 'pseudomallei group', 'Teredinibacter turnerae', 'Myxococcus', 'Plenodomus lingam/Leptosphaeria maculans species complex', 'Purpureocillium', 'Chloroflexales', 'Saccharomycetales', 'Neonectria', 'Ecdysozoa', 'Hypocreaceae', 'Orbiliaceae', 'Gelatoporiaceae', 'Nectriaceae', 'Eleftheria', 'Gelatoporia', 'Candidatus Entotheonellaceae', 'Sandaracinus', 'Orbilia oligospora', 'Myxococcaceae', 'Bacilli', 'Cyanodermella', 'Shigella', 'Macrophomina', 'Metarhizium', 'Brevibacillus', 'Actinomadura', 'Cellvibrionaceae', 'Labrenzia', 'Basidiomycota', 'Kamptonema', 'Polyporaceae', 'Thermobifida', 'Nocardia', 'Xanthomonadales', 'Plenodomus', 'Mycobacterium avium complex (MAC)', 'Herbaspirillum', 'Myxococcia', 'Pleosporales', 'Neoptera', 'Sphaerotilaceae', 'Lysobacteraceae', 'Saccharopolyspora', 'Listeria', 'Microascaceae', 'Tapinella', 'Alternaria sect. Alternaria', 'Beauveria', 'Acinetobacter', 'Trinickia', 'Burkholderiales', 'Amycolatopsis japonica group', 'Streptantibioticus', 'Oscillatoriophycideae', 'Strobilurus', 'Candidatus Entotheonellales', 'Dothideales', 'Geminicoccales', 'Streptosporangiaceae', 'Pseudomonadota', 'Microcystaceae', 'Phialomyces', 'Photobacterium', 'Dipodascales', 'Aspergillus subgen. Aspergillus', 'Ascomycota', 'Actinomycetota', 'Mycetohabitans', 'Kiritimatiellota', 'environmental samples', 'Bacillota', 'Massarineae', 'Halomonadaceae', 'Saccharomycotina', 'Lactobacillaceae', 'Bionectriaceae', 'Actinoplanes', 'Bipolaris', 'Streptococcus', 'Comamonadaceae', 'Klebsiella/Raoultella group', 'Cyanobacteria', 'Didymellaceae', 'Chloroflexi bacterium', 'Lactobacillales', 'Opitutia', 'Nocardiopsidaceae', 'Pseudonocardiaceae', 'Aureobasidium', 'other sequences'})

    df = process_list(names)
    # Ensure column order
    cols = ["supplied_name","matched_name","tax_id","rank","kingdom","phylum","class","order","family","genus","species","notes"]
    df = df[cols]
    df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(df)} rows to {args.output}")


if __name__ == "__main__":
    main()