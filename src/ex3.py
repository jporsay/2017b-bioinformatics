import argparse
import os
from collections import namedtuple

from Bio import AlignIO

BlastConfig = namedtuple("BlastConfig", "cmd, db, evalue")


class Params(object):
    def __init__(self):
        self.source = None

    def validate(self):
        if not os.path.isfile(self.source):
            raise FileNotFoundError("%s doesn't exist" % (self.source,))


def parse_args() -> Params:
    parser = argparse.ArgumentParser(description='Performs MSA with given FASTA sequences')
    parser.add_argument("--source", help="Fasta source file", required=True)
    return parser.parse_args(namespace=Params())


def main():
    args = parse_args()
    args.validate()
    alignments = AlignIO.read(args.source, "fasta")
    print(alignments)


if __name__ == '__main__':
    main()
