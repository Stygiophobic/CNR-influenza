#singularity shell /srv/nfs/ngs-stockage/NGS_Virologie/NEXTSTRAIN/nextstrainV3.simg
#cp /srv/nfs/ngs-stockage/NGS_Virologie/hcl-vir-ngs/CNRVI/2019_2020/gisaid_epiflu_uploader_v113_surveillance20190109.xls ~/git/CNR-influenza/data/last_gisaid_xls.xls

dataset = ['H1N1S4','H1N1S6','H3N2S4','H3N2S6','BS4','BS6']

rule all:
    input:
        H1N1_S4 = "temp_data/H1N1S4.tsv",
        H1N1_S6 = "temp_data/H1N1S6.tsv",
        H3N2_S4 = "temp_data/H3N2S4.tsv",
        H3N2_S6 = "temp_data/H3N2S6.tsv",
        B_S4 = "temp_data/BS4.tsv",
        B_S6 = "temp_data/BS6.tsv",
        B_S4t = "temp_data/BS4.nwk",
        B_S6t = "temp_data/BS6.nwk" ,  
        H1N1_S4t = "temp_data/H1N1S4.nwk",
        H1N1_S6t = "temp_data/H1N1S6.nwk",
        H3N2_S4t = "temp_data/H3N2S4.nwk",
        H3N2_S6t = "temp_data/H3N2S6.nwk",
        auspice_tree_H1N1_S4 = "auspice/CNR-influenza_H1N1S4.json",
        auspice_tree_H1N1_S6 = "auspice/CNR-influenza_H1N1S6.json",
        auspice_tree_H3N2_S4 = "auspice/CNR-influenza_H3N2S4.json",
        auspice_tree_H3N2_S6 = "auspice/CNR-influenza_H3N2S6.json",
        auspice_tree_B_S4 = "auspice/CNR-influenza_BS4.json",
        auspice_tree_B_S6 = "auspice/CNR-influenza_BS6.json",
        node_dataH14 = "temp_data/H1N1S4_nt_muts.json",
        node_dataH16 = "temp_data/H1N1S6_nt_muts.json",
        node_dataH34 = "temp_data/H3N2S4_nt_muts.json",
        node_dataH36 = "temp_data/H3N2S6_nt_muts.json",
        node_dataB4 = "temp_data/BS4_nt_muts.json",
        node_dataB6 = "temp_data/BS6_nt_muts.json",
        node_dataH14a = "temp_data/H1N1S4_aa_muts.json",
        node_dataH16a = "temp_data/H1N1S6_aa_muts.json",
        node_dataH34a = "temp_data/H3N2S4_aa_muts.json",
        node_dataH36a = "temp_data/H3N2S6_aa_muts.json",
        node_dataB4a = "temp_data/BS4_aa_muts.json",
        node_dataB6a = "temp_data/BS6_aa_muts.json"


rule get_last_data:
    input:
    output:

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
        H1N1_S4 = "temp_data/H1N1S4.tsv",
        H1N1_S6 = "temp_data/H1N1S6.tsv",
        H3N2_S4 = "temp_data/H3N2S4.tsv",
        H3N2_S6 = "temp_data/H3N2S6.tsv",
        B_S4 = "temp_data/BS4.tsv",
        B_S6 = "temp_data/BS6.tsv"
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
        filter_fasta = rules.augur_filter.output.filtered_seq ,
        ref_seq = "config/{subset}.fasta"
    output:
        align_fasta = "temp_data/{subset}_align.fasta"
    shell:
        "augur align "
        "--sequences {input.filter_fasta} "
        "--reference-sequence {input.ref_seq} "
        "--remove-reference "
        "--fill-gaps "
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

rule augur_ancestral:
    input:
        tree = rules.augur_refine.output.tree,
        alignment = rules.augur_align.output.align_fasta
    output:
        node_data = "temp_data/{subset}_nt_muts.json"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data}
        """

rule augur_translate:
    input:
        tree = rules.augur_refine.output.tree,
        node_data = rules.augur_ancestral.output.node_data,
        reference = "config/{subset}.gb"
    output:
        node_data = "temp_data/{subset}_aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """

rule export:
    input:
        tree = rules.augur_refine.output.tree,
        metadata = "temp_data/{subset}.tsv",
        nt_muts = rules.augur_ancestral.output.node_data,
        aa_muts = rules.augur_translate.output.node_data,
        branch_lengths = rules.augur_refine.output.node_data
        #auspice_config = "config/auspice_config.json"
    output:
        auspice_json = "auspice/CNR-influenza_{subset}.json",
    shell:
        "augur export v2 "
        "--tree {input.tree} "
        "--metadata {input.metadata} "
        "--node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} "
        #"--auspice-config {input.auspice_config} "
        "--output {output.auspice_json} "
