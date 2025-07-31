from Bio import SeqIO
import argparse


def genbank_to_bed(genbank_file, bed_file):
    with open(bed_file, "w") as bed:
        for record in SeqIO.parse(genbank_file, "genbank"):
            for feature in record.features:
                record.id = record.name
                if feature.type != "source" and feature.type != "primer_bind":
                    for quali in feature.qualifiers:
                        if quali == 'label':
                            start = feature.location.start
                            end = feature.location.end
                            name = feature.qualifiers['label'][0]
                            if feature.location.strand < 0:
                              strand = "-"
                            else:
                             strand = "+"
                            if (end-start) > 50:
                                bed.write(f"{record.id}\t{start}\t{end}\t{name}\t0\t{strand}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GenBank to BED")
    parser.add_argument("genbank_file", help="Input GenBank file")
    parser.add_argument("bed_file", help="Output BED file")
    args = parser.parse_args()

    genbank_to_bed(args.genbank_file, args.bed_file)