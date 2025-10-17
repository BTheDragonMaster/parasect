from parasect.core.parsing import parse_fasta_file
from sys import argv
from sqlalchemy import create_engine
from sqlalchemy.orm import Session
from parasect.database.query_database import get_proteins_from_synonym
if __name__ == "__main__":
    engine = create_engine(f"sqlite:///{argv[1]}")
    fasta = parse_fasta_file(argv[2])
    bacteria = 0
    fungi = 0
    other = 0
    unknown = 0
    with Session(engine) as session:
        for seq_id in fasta:
            prot_id = seq_id.split('\t')[2]
            protein = get_proteins_from_synonym(session, prot_id)
            if protein:
                protein = protein[0]
                if protein.taxonomy.domain == "Bacteria":
                    bacteria += 1
                elif protein.taxonomy.kingdom == "Fungi":
                    fungi += 1
                else:
                    other += 1
            else:
                unknown += 1

    print(bacteria, fungi, other, unknown)


