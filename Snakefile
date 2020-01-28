#singularity shell ../singularity_nextstrain/nextstrainV3.simg

rule all:
    input:
        xls_file = "data/last_gisaid_xls.xls",
        H1N1_S4 = "temp_data/H1N1_S4.tsv",
        H1N1_S6 = "temp_data/H1N1_S6.tsv",
        H3N2_S4 = "temp_data/H3N2_S4.tsv",
        H3N2_S6 = "temp_data/H3N2_S6.tsv",
        B_S4 = "temp_data/B_S4.tsv",
        B_S6 = "temp_data/B_S6.tsv",
        H1N1_S4a = "temp_data/H1N1_S4_align.fasta",
        H1N1_S6a = "temp_data/H1N1_S6_align.fasta",
        H3N2_S4a = "temp_data/H3N2_S4_align.fasta",
        H3N2_S6a = "temp_data/H3N2_S6_align.fasta",
        B_S4a = "temp_data/B_S4_align.fasta",
        B_S6a = "temp_data/B_S6_align.fasta"                  

rule xls_to_fasta_csv:
    input:
        xls_file = rules.all.input.xls_file
    output:
        metadata_raw = "temp_data/metadata_raw.csv",
        fasta_seq = "temp_data/sequences.fasta"      
    shell:
        "script/process_xls.py {input} {output.fasta_seq} {output.metadata_raw}"   

    
rule make_metadata:
    input:
        csv_file = rules.xls_to_fasta_csv.output.metadata_raw
        #csv_file = "temp_data/metadata_raw.csv"
    output:
        #metadata = "temp_data/{subset}.tsv"
        H1N1_S4 = "temp_data/H1N1_S4.tsv",
        H1N1_S6 = "temp_data/H1N1_S6.tsv",
        H3N2_S4 = "temp_data/H3N2_S4.tsv",
        H3N2_S6 = "temp_data/H3N2_S6.tsv",
        B_S4 = "temp_data/B_S4.tsv",
        B_S6 = "temp_data/B_S6.tsv"
    shell:
        "Rscript script/make_metadata.R {input} "
        
rule augur_filter:
    input:
        seq_file = rules.xls_to_fasta_csv.output.fasta_seq ,
        meta_subtype = "temp_data/{subset}.tsv"

    output:
        filtered_seq = "temp_data/{subset}.fasta"  
    shell:
        "augur filter  "
        "--sequences {input.seq_file} "
        "--metadata {input.meta_subtype}  "
        "--output {output.filtered_seq} " 

rule augur_align:
    input:
        filter_fasta = rules.augur_filter.output.filtered_seq
    output:
        align_fasta = "temp_data/{subset}_align.fasta"
    shell:
        "augur align "
        "--sequences {input} "
        "--output {output} "
            


