import argparse
import os
from collections import namedtuple

from Bio import AlignIO
from Bio.Blast import NCBIWWW, NCBIXML

BlastConfig = namedtuple("BlastConfig", "cmd, db, evalue")


class Params(object):
    def __init__(self):
        self.source = None
        self.dest = None

    def validate(self):
        if not os.path.isfile(self.source):
            raise FileNotFoundError("%s doesn't exist" % (self.source,))


def parse_args() -> Params:
    parser = argparse.ArgumentParser(description='Performs MSA with the best 10 matches from a BLAST query for a given'
                                                 'FASTA sequence')
    parser.add_argument("--source", help="Fasta source file", required=True)
    parser.add_argument("--dest", help="Output file", required=True)
    return parser.parse_args(namespace=Params())


def main():
    args = parse_args()
    args.validate()
    config = BlastConfig("blastx", "swissprot", 0.0001)
    with open(args.source) as fasta_input:
        handle = NCBIWWW.qblast(config.cmd, config.db, fasta_input.read(), expect=config.evalue)
        top_10_alignments = next(NCBIXML.parse(handle)).alignments[:10]
        original_alignment = AlignIO.read(args.source, "fasta")
        print(original_alignment)
        for record in top_10_alignments:
            print(record)


if __name__ == '__main__':
    main()
