import argparse
import os

from Bio import SeqIO


class Params(object):
    def __init__(self):
        self.source = None
        self.dest = None

    def validate(self):
        if not os.path.isfile(self.source):
            raise FileNotFoundError("%s doesn't exist" % (self.source,))

    def __str__(self):
        return "source=%s dest=%s" % (self.source, self.dest)


def parse_args() -> Params:
    parser = argparse.ArgumentParser(description='Converts genebank sequences to fasta')
    parser.add_argument("--source", help="Genebank source file", required=True)
    parser.add_argument("--dest", help="Fasta target file", required=True)
    return parser.parse_args(namespace=Params())


def main():
    args = parse_args()
    args.validate()
    count = SeqIO.convert(args.source, "genbank", args.dest, "fasta")
    print("Converted %i records" % count)


if __name__ == '__main__':
    main()
