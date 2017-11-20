import argparse
import os
import shutil
from io import StringIO

from blast_utils import blast_local, blast_remote, BlastConfig


class Params(object):
    def __init__(self):
        self.source = None
        self.dest = None
        self.blast = None

    def validate(self):
        if not os.path.isfile(self.source):
            raise FileNotFoundError("%s doesn't exist" % (self.source,))


def parse_args() -> Params:
    parser = argparse.ArgumentParser(description='Runs a BLAST query using data from a local FASTA file')
    parser.add_argument("--source", help="Fasta source file", required=True)
    parser.add_argument("--dest", help="Output file", required=True)
    parser.add_argument("--blast", help="Use remote or local blast", required=True, choices=["local", "remote"])
    return parser.parse_args(namespace=Params())


def process_results(args, handle: StringIO):
    with open(args.dest, "w") as out:
        handle.seek(0)
        shutil.copyfileobj(handle, out)
        # for blast_record in NCBIXML.parse(handle):
        #     out.write("Alignments: {}\n".format(len(blast_record.alignments)))
        #     for alignment in blast_record.alignments:
        #         out.write(f"\n\n>>>>> Alignment {alignment.title}\nId: {alignment.hit_id}\n")
        #         for hsp in alignment.hsps:
        #             out.write(f'Score {hsp.score} ({hsp.bits} bits) expectation {hsp.expect}, '
        #                       f'alignment length {hsp.align_length}\n')
        #             out.write(f'Query\t{hsp.query_start}\t{hsp.match} {hsp.query_end}\n')
        #             out.write(f'Subj\t{hsp.sbjct_start}\t\t{hsp.sbjct} {hsp.sbjct_end}\n')
        #             out.write('\n')


def main():
    args = parse_args()
    args.validate()
    config = BlastConfig("blastx", "swissprot", 0.0001)
    blast_fn = blast_local if "local" in args.blast else blast_remote
    output = blast_fn(args.source, config)
    process_results(args, output)
    print("Results can be found inside {}".format(args.dest))


if __name__ == '__main__':
    main()
