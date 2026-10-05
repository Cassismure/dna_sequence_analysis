from Bio import Entrez, SeqIO
import csv
import time

Entrez.email = "okeleyetomilola@gmail.com"

TARGETS = [
    ("Homo sapiens", "TP53"),
    ("Homo sapiens", "BRCA1"),
    ("Escherichia coli", "lacZ"),
    ("Saccharomyces cerevisiae", "ACT1"),
    ("Drosophila melanogaster", "white"),
]


def fetch_gene_sequence(organism: str, gene_name: str):
    """
    Search NCBI's nucleotide database for an mRNA record matching
    this organism + gene name, and fetch the first hit as FASTA.
    Returns a Biopython SeqRecord, or None if nothing was found.
    """
    query = f"{organism}[Organism] AND {gene_name}[Gene Name] AND biomol_mrna[PROP]"
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=1)
    search_result = Entrez.read(handle)
    handle.close()

    id_list = search_result["IdList"]
    if not id_list:
        print(f"  No result for {organism} / {gene_name}")
        return None

    fetch_handle = Entrez.efetch(db="nucleotide", id=id_list[0], rettype="fasta", retmode="text")
    record = SeqIO.read(fetch_handle, "fasta")
    fetch_handle.close()
    return record


def main():
    rows = []
    for i, (organism, gene) in enumerate(TARGETS, start=1):
        print(f"Fetching {organism} / {gene} ...")
        record = fetch_gene_sequence(organism, gene)
        if record is not None:
            # Some mRNA records are long — trim to first 300bp so the
            # dataset stays comparable in scale to the synthetic version
            seq = str(record.seq)[:300]
            rows.append((f"REAL{i:03d}", organism, seq))
        time.sleep(0.4) 

    with open("data/sequences.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sequence_id", "organism", "sequence"])
        writer.writerows(rows)

    print(f"\nSaved {len(rows)} real sequences to data/sequences.csv")


if __name__ == "__main__":
    main()
