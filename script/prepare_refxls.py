#!/usr/bin/env python3

#singularity shell /srv/nfs/ngs-stockage/NGS_Virologie/NEXTSTRAIN/nextstrainV3.simg

import pandas as pd
import sys

#xls_file="/srv/nfs/ngs-stockage/NGS_Virologie/hcl-vir-ngs/CNRVI/2019_2020/gisaid_epiflu_isolates_H1N1_20200203.xlsx"
#csv_file="/srv/nfs/ngs-stockage/NGS_Virologie/NEXTSTRAIN/data.csv"
xls_file=sys.argv[1]
#csv_file=sys.argv[2]

meta_xls = pd.read_excel(xls_file)
#meta_xls.to_csv (csv_file, index = None,sep=";" ,header=True)
meta_xls.to_csv ("temp.csv", index = None,sep=";" ,header=True)