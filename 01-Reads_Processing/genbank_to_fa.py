from Bio import SeqIO
import argparse

def genbank_to_fasta(genbank_file, fasta_file):
    with open(genbank_file) as input_handle:
        with open(fasta_file, "w") as output_handle:
            for record in SeqIO.parse(input_handle, "genbank"):
                record.id = record.name
                SeqIO.write(record, output_handle, "fasta-2line")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GenBank to FASTA")
    parser.add_argument("genbank_file", help="GenBank file")
    parser.add_argument("fasta_file", help="Output FASTA file")
    args = parser.parse_args()

    genbank_to_fasta(args.genbank_file, args.fasta_file)