#!/bin/sh
set -x

getorf -sequence tsc1.gb -outseq ng_012386.orf
echo "Wrote ORFs to ng_012386.orf"
