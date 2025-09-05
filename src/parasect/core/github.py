from typing import Any
import os
from sys import argv
import json
import requests
from dotenv import load_dotenv
from dataclasses import dataclass
from enum import Enum

from sqlalchemy.orm import Session

from parasect.database.query_database import get_substrates_from_name
from parasect.core.chem import is_same_molecule


load_dotenv()
GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN")
REPO = "BTheDragonMaster/parasect"
API_URL = f"https://api.github.com/repos/{REPO}/issues"


class EntryType(Enum):
    NEW = 1  # New domain entry, not seen in dataset
    CORRECTION = 2  # Update to domain in the dataset
    DUPLICATE = 3  # Domain entered with the SAME name but a DIFFERENT sequence

    @staticmethod
    def from_string(string: str) -> "EntryType":
        string_to_enum = {"new_entry": EntryType.NEW,
                          "correction": EntryType.CORRECTION,
                          "duplicate_entry": EntryType.DUPLICATE}
        if string in string_to_enum:
            return string_to_enum[string]
        else:
            raise ValueError(f"Unknown entry type: {string}. \
            Must be one of 'new_entry', 'correction', or 'duplicate_entry'")

    @staticmethod
    def to_header(entry_type: "EntryType") -> str:
        enum_to_header_string = {EntryType.NEW: "New domain entries",
                                 EntryType.CORRECTION: "Domain corrections",
                                 EntryType.DUPLICATE: "Domain sequence revision"}

        return enum_to_header_string[entry_type]


@dataclass
class ProteinEntry:
    name: str
    sequence: str


@dataclass
class SubstrateEntry:
    name: str
    smiles: str
    new: bool

    def __eq__(self, other):
        if type(self) == type(other) and self.name == other.name and is_same_molecule(self.smiles, other.smiles):
            return True

        return False

    def __hash__(self):
        return hash(self.name)


@dataclass
class DomainEntry:
    name: str
    sequence: str
    signature: str
    extended_signature: str
    substrates: list[SubstrateEntry]
    entry_type: EntryType


def submit_github_issue_protein(
    protein_entry, 
    domains, 
    new_substrates, 
    annotation_type: EntryType,
    orcid: str = "",
    references: list[dict[str, str]] = []
):
    url = f"https://api.github.com/repos/{REPO}/issues"
    headers = {"Authorization": f"token {GITHUB_TOKEN}"}

    title = f"[{EntryType.to_header(annotation_type)}] {protein_entry.name}"
    body_dict = {
        "orcid": orcid,
        "references": references,
        "protein_name": protein_entry.name,
        "protein_sequence": protein_entry.sequence,
        "new_substrates": [{"name": substrate.name, "smiles": substrate.smiles} for substrate in new_substrates],
        "domains": [{"name": domain.name,
                     "sequence": domain.sequence,
                     "signature": domain.signature,
                     "extended_signature": domain.extended_signature,
                     "substrates": [s.name for s in domain.substrates]} for domain in domains]
    }

    body = json.dumps(body_dict, indent=2)

    issue = {"title": title, "body": body}
    response = requests.post(url, json=issue, headers=headers)

    # evaluate response code and throw error if not good
    if response.status_code != 201:
        raise ValueError(f"Failed to create GitHub issue: {response.text}")

    return response.json()


def submit_github_issues(
    session: Session, 
    protein_to_entry: dict[str, dict[str, Any]],
    orcid: str,
    references: list[dict[str, str]]
) -> None:
    """Submit GitHub new annotations, one per protein per annotation type, to BTheDragonMaster/parasect as issues

    :param session: Database session
    :type session: Session
    :param protein_to_entry: dictionary of protein names as keys and protein and domain annotations as (nested) values
    :type protein_to_entry: dict[str, dict[str, Any]]
    :param orcid: ORCID ID of the user submitting the annotations
    :type orcid: str
    :param references: list of references to include in the issue
    :type references: list[dict[str, str]]
    """
    for protein_name, protein_info in protein_to_entry.items():
        protein_entry = ProteinEntry(protein_name, protein_info["sequence"])

        corrections: list[DomainEntry] = []
        new_entries: list[DomainEntry] = []
        duplicate_entries: list[DomainEntry] = []

        substrates_corrections: set[SubstrateEntry] = set()
        substrates_new_entries: set[SubstrateEntry] = set()
        substrates_duplicates: set[SubstrateEntry] = set()

        domain_to_entry = protein_info["domains"]

        for domain, domain_info in domain_to_entry.items():

            annotation_type = EntryType.from_string(domain_info["annotation_type"])
            domain_substrates: list[SubstrateEntry] = []
            substrates = domain_info["substrates"]
            for substrate in substrates:
                print(substrate)

                if not get_substrates_from_name(session, substrate["name"]):
                    new = True
                else:
                    new = False

                substrate_entry = SubstrateEntry(substrate["name"],
                                                 substrate["smiles"],
                                                 new)

                domain_substrates.append(substrate_entry)

            domain_entry = DomainEntry(domain,
                                       domain_info["sequence"],
                                       domain_info["signature"],
                                       domain_info["extended_signature"],
                                       domain_substrates, annotation_type)
            if domain_entry.entry_type == EntryType.NEW:

                new_entries.append(domain_entry)
                for substrate in domain_substrates:
                    if substrate.new:
                        substrates_new_entries.add(substrate)
            elif domain_entry.entry_type == EntryType.CORRECTION:
                corrections.append(domain_entry)
                for substrate in domain_substrates:
                    if substrate.new:
                        substrates_corrections.add(substrate)
            elif domain_entry.entry_type == EntryType.DUPLICATE:
                duplicate_entries.append(domain_entry)
                for substrate in domain_substrates:
                    if substrate.new:
                        substrates_duplicates.add(substrate)
            else:
                raise ValueError(f"Unrecognised entry type: {domain_entry.entry_type.name}")

        if corrections:
            submit_github_issue_protein(protein_entry, corrections, substrates_corrections, EntryType.CORRECTION, orcid, references)
        if new_entries:
            submit_github_issue_protein(protein_entry, new_entries, substrates_new_entries, EntryType.NEW, orcid, references)
        if duplicate_entries:
            submit_github_issue_protein(protein_entry, duplicate_entries, substrates_duplicates, EntryType.DUPLICATE, orcid, references)


def fetch_github_issues(token, filter_strings=None, state="all"):
    """
    Fetches all issues from the repo and filters by title content.

    :param token: GitHub personal access token
    :param filter_strings: list of substrings to match in the title
    :param state: 'open', 'closed', or 'all'
    :return: list of matching issue dicts
    """
    headers = {"Authorization": f"token {token}"}
    page = 1
    all_matches = []

    while True:
        params = {"state": state, "per_page": 100, "page": page}
        response = requests.get(API_URL, headers=headers, params=params)
        if response.status_code != 200:
            raise RuntimeError(f"Error fetching issues: {response.status_code}, {response.text}")

        issues = response.json()
        if not issues:
            break  # No more pages

        if filter_strings:
            for issue in issues:
                title = issue["title"]
                if any(f in title for f in filter_strings):
                    all_matches.append(issue)
        else:
            all_matches.extend(issues)

        page += 1

    return all_matches


def issues_to_files(out_dir: str) -> None:
    """Write GitHub issues to folder structure for database processing

    :param out_dir: output directory
    :type out_dir: str
    """
    keywords = ["New domain entries", "Domain corrections", "Domain sequence revision"]
    issues = fetch_github_issues(GITHUB_TOKEN, filter_strings=keywords)

    if issues:
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        for issue in issues:
            folder_title = f"issue_{issue['number']}"
            issue_folder = os.path.join(out_dir, folder_title)
            if not os.path.exists(issue_folder):
                os.mkdir(issue_folder)

            data = json.loads(issue["body"])
            if data["new_substrates"]:
                smiles_path = os.path.join(issue_folder, 'smiles.tsv')
                with open(smiles_path, 'w') as smiles_out:
                    smiles_out.write("substrate\tsmiles\n")
                    for substrate in data["new_substrates"]:
                        smiles_out.write(f"{substrate['name']}\t{substrate['smiles']}\n")
            protein_path = os.path.join(issue_folder, "proteins.fasta")
            with open(protein_path, 'w') as protein_out:
                protein_out.write(f">{data['protein_name']}\n{data['protein_sequence']}\n")

            domain_path = os.path.join(issue_folder, "domains.fasta")
            signature_path = os.path.join(issue_folder, "signatures.fasta")
            extended_path = os.path.join(issue_folder, "extended_signatures.fasta")
            parasect_path = os.path.join(issue_folder, "parasect_dataset.txt")
            with open(domain_path, 'w') as domains_out, open(parasect_path, 'w') as parasect_out, \
                    open(signature_path, 'w') as signatures_out, open(extended_path, 'w') as extended_out:
                parasect_out.write("domain_id\tsequence\tspecificity\n")
                for domain in data["domains"]:
                    domains_out.write(f">{domain['name']}\n{domain['sequence']}\n")
                    signatures_out.write(f">{domain['name']}\n{domain['signature']}\n")
                    extended_out.write(f">{domain['name']}\n{domain['extended_signature']}\n")
                    parasect_out.write(f"{domain['name']}\t{domain['sequence']}\t{'|'.join(domain['substrates'])}\n")


if __name__ == "__main__":
    issues_to_files(argv[1])



