from collections import namedtuple
from io import StringIO
from tempfile import NamedTemporaryFile

import os
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline

BlastConfig = namedtuple("BlastConfig", "cmd, db, evalue")


def blast_remote(source: str, config: BlastConfig) -> StringIO:
    print("Running with remote BLAST")
    with open(source) as fasta_input:
        handle = NCBIWWW.qblast(config.cmd, config.db, fasta_input.read(), expect=config.evalue)
    return handle


def blast_local(source: str, config: BlastConfig) -> StringIO:
    print("Running with local BLAST")
    ntf = NamedTemporaryFile(delete=False)
    temp_name = ntf.name
    ntf.close()
    cli = NcbiblastxCommandline(cmd=config.cmd, db=config.db, evalue=config.evalue, query=source, out=temp_name, outfmt=5)
    cli()
    out = StringIO()
    with open(temp_name, "r") as handle:
        out.write(handle.read())
    os.remove(temp_name)
    return out
