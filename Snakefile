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

rule make_metadata:
    input:
        csv_file = rules.xls_to_fasta_csv.output.metadata_raw
    output:
        metadata = "temp_data/{subset}.tsv"
        #H1N1_S4 = "temp_data/H1N1_S4.tsv",
        #H1N1_S6 = "temp_data/H1N1_S6.tsv",
        #H3N2_S4 = "temp_data/H3N2_S4.tsv",
        #H3N2_S6 = "temp_data/H3N2_S6.tsv",
        #B_S4 = "temp_data/B_S4.tsv",
        #B_S6 = "temp_data/B_S6.tsv"
    message:
        "production des fichiers de metadata"
    script:
        "script/make_metadata.R"
        #"Rscript script/make_metadata.R {input.csv_file}" 
        #"output.H1N1_S4.tsv H1N1_S6.tsv"
        #"H3N2_S4.tsv H3N2_S6.tsv"
        #"B_S4.tsv B_S6.tsv"    
        #"{output.H1N1_S4} {output.H1N1_S6k}"
        #"{output.H3N2_S4} {output.H3N2_S6}"
        #"{output.B_S4} {output.B_S6}"
                  
rule augur_filter:
    input:
        seq_file = rules.xls_to_fasta_csv.output.fasta_seq ,
        meta_subtype = "{subset}.tsv"
    output:
        filtered_seq = "{subset}.fasta"
    message:
        "MIAM MIAM"    
    shell:
        "augur filter \ "
        "--sequence {input.seq_file} \ "
        "--metadata {input.meta_subtype} \ "
        "--output {output.filtered_seq} "
