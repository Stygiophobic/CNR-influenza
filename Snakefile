#!/usr/bin/env python3

#singularity shell ../singularity_nextstrain/nextstrainV3.simg

rule xls_to_fasta_csv:
    input:
        xls_file = "data/last_gisaid_xls.xls"
    output:
        metadata_raw = "temp_data/metadata_raw.csv",
        fasta_seq = "temp_data/sequences.fasta"
    message:
        "ALLO MDR KEK"        
    shell:
        "script/process_xls.py {input} {output.fasta_seq} {output.metadata_raw}"        


