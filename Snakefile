#singularity shell ../singularity_nextstrain/nextstrainV3.simg

rule all:
    input:
        auspice_tree_H1N1_S4 = "auspice/CNR-influenza_H1N1_S4_tree.json",
        auspice_meta_H1N1_S4 = "auspice/CNR-influenza_H1N1_S4_meta.json",
        auspice_tree_H1N1_S6 = "auspice/CNR-influenza_H1N1_S6_tree.json",
        auspice_meta_H1N1_S6 = "auspice/CNR-influenza_H1N1_S6_meta.json",
        auspice_tree_H3N2_S4 = "auspice/CNR-influenza_H3N2_S4_tree.json",
        auspice_meta_H3N2_S4 = "auspice/CNR-influenza_H3N2_S4_meta.json",
        auspice_tree_H3N2_S6 = "auspice/CNR-influenza_H3N2_S6_tree.json",
        auspice_meta_H3N2_S6 = "auspice/CNR-influenza_H3N2_S6_meta.json",
        auspice_tree_B_S4 = "auspice/CNR-influenza_B_S4_tree.json",
        auspice_meta_B_S4 = "auspice/CNR-influenza_B_S4_meta.json",
        auspice_tree_B_S6 = "auspice/CNR-influenza_B_S6_tree.json",
        auspice_meta_B_S6 = "auspice/CNR-influenza_B_S6_meta.json"     


rule xls_to_fasta_csv:
    input:
        xls_file = "data/last_gisaid_xls.xls"
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
            
rule augur_raw_tree:
    input:
        align_data = rules.augur_align.output.align_fasta
    output:
        raw_tree = "temp_data/{subset}_raw.nwk"
    shell:
        "augur tree "
        "--alignment {input} "
        "--output {output} "
                

rule augur_refine:
    input:
        tree = rules.augur_raw_tree.output.raw_tree,
        alignment = rules.augur_align.output.align_fasta,
        meta  = "temp_data/{subset}.tsv"
    output:
        tree = "temp_data/{subset}.nwk",
        node_data = "temp_data/{subset}_branch_lengths.json"
    shell:
        "augur refine "
        "--tree {input.tree} "
        "--alignment {input.alignment} "
        "--metadata {input.meta} "
        "--timetree "
        "--output-tree {output.tree} "
        "--output-node-data {output.node_data} "

rule augur_export:
    input:
        tree = rules.augur_refine.output.tree,
        meta  = "temp_data/{subset}.tsv",
        branch_lengths = rules.augur_refine.output.node_data,
        auspice_config = "config/auspice_config.json"
    output:
        auspice_tree = "auspice/CNR-influenza_{subset}_tree.json",
        auspice_meta = "auspice/CNR-influenza_{subset}_meta.json"
    shell:
        "augur export v1 "
        "--tree {input.tree} "
        "--metadata {input.meta} "
        "--node-data {input.branch_lengths} "
        "--auspice-config {input.auspice_config} "
        "--output-tree {output.auspice_tree} "
        "--output-meta {output.auspice_meta} "
  