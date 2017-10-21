import argparse
import os
from collections import namedtuple
from tempfile import NamedTemporaryFile

from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline

BlastConfig = namedtuple("BlastConfig", "cmd, db, evalue")


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


def blast_remote(args: Params, config: BlastConfig):
    print("Running with remote BLAST")
    with open(args.source) as fasta_input:
        handle = NCBIWWW.qblast(config.cmd, config.db, fasta_input.read(), expect=config.evalue)
        process_results(args, handle)


def blast_local(args: Params, config: BlastConfig):
    print("Running with local BLAST")
    ntf = NamedTemporaryFile(delete=False)
    temp_name = ntf.name
    ntf.close()
    cli = NcbiblastxCommandline(cmd=config.cmd, db=config.db, evalue=config.evalue, query=args.source, out=temp_name, outfmt=5)
    cli()
    with open(temp_name, "r") as handle:
        process_results(args, handle)
    os.remove(temp_name)


def process_results(args, handle):
    with open(args.dest, "w") as out:
        for blast_record in NCBIXML.parse(handle):
            out.write("Alignments: {}\n".format(len(blast_record.alignments)))
            for alignment in blast_record.alignments:
                out.write("\n\n>>>>> Alignment {}\n".format(alignment.title))
                for hsp in alignment.hsps:
                    out.write(str(hsp))
                    out.write("\n")


def main():
    args = parse_args()
    args.validate()
    config = BlastConfig("blastx", "swissprot", 0.0001)
    blast_fn = blast_local if "local" in args.blast else blast_remote
    blast_fn(args, config)
    print("Results can be found inside {}".format(args.dest))


if __name__ == '__main__':
    main()
