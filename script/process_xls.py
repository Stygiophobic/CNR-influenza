#!/usr/bin/env python3

import pandas as pd
import sys

xls_file=sys.argv[1]
seq_file=sys.argv[2]
csv_file=sys.argv[3]

#extraire la feuille "Isolates" et la convertir en .csv
meta_xls = pd.read_excel(xls_file, sheet_name='Isolates')
meta_xls.to_csv (csv_file, index = None,sep=";" ,header=True)

#extraire les sequences depuis la feuille xls "Sequences"
sequences_xls = pd.read_excel(xls_file, sheet_name='Sequences',header=None)
nrow=len(sequences_xls)
fasta_file=open(seq_file,"w")
for i in range(0,nrow):
    fasta_file.write(sequences_xls.iat[i, 0]+"\n")
fasta_file.close()    