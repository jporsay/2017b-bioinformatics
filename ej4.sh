#!/bin/sh
set -x

./venv/bin/python ./src/ex4.py --source ./ej2_out.txt --pattern "Cytosolic carboxypeptidase" --savematch match
