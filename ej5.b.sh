#!/bin/sh
set -x

patmatmotifs -full ./tsc1.fasta -outfile tsc1.patmatmotifs
echo "Result can be found inside tsc1.patmatmotifs"
