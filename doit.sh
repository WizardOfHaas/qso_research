#!/bin/bash
rm plots/* -rf
python new_abs.py possible_bal.csv
pdfunite plots/*.pdf plots/all.pdf