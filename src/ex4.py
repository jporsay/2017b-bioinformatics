import argparse
import os

from Bio.Blast import NCBIXML
from Bio import Entrez


class Params(object):
    def __init__(self):
        self.source = None
        self.pattern = None
        self.savematch = None

    def validate(self):
        if not os.path.isfile(self.source):
            raise FileNotFoundError("%s doesn't exist" % (self.source,))


def parse_args() -> Params:
    parser = argparse.ArgumentParser(description='Runs a BLAST query using data from a local FASTA file')
    parser.add_argument("--source", help="Blast output file", required=True)
    parser.add_argument("--pattern", help="Pattern to search for", required=True)
    parser.add_argument("--savematch", help="Save match fasta file", required=True)
    return parser.parse_args(namespace=Params())


def process_blast(source, dest, pattern):
    match = False
    hit = 0
    with open(source, "r") as blast_handle:
        for blast_record in NCBIXML.parse(blast_handle):
            for alignment in blast_record.alignments:
                if pattern in alignment.title:
                    match = True
                    Entrez.email = "jorsay@itba.edu.ar"
                    handle = Entrez.efetch(db="protein", id=alignment.accession, rettype="fasta", retmode="text")
                    dest_file = dest + "-" + str(hit) + ".fasta"
                    with open(dest_file, "w") as out:
                        out.write(handle.read())
                    print("Found a match. Saved it in", dest_file)
                    hit = hit + 1
    if not match:
        print("Couldn't find a match :(")


def main():
    args = parse_args()
    args.validate()
    # process_blast(args.source, args.pattern)
    process_blast(args.source, args.savematch, args.pattern)
    print("Results can be found inside {}".format(args.savematch))


if __name__ == '__main__':
    main()
