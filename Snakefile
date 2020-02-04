#singularity shell /srv/nfs/ngs-stockage/NGS_Virologie/NEXTSTRAIN/nextstrainV3.simg
#cp /srv/nfs/ngs-stockage/NGS_Virologie/hcl-vir-ngs/CNRVI/2019_2020/gisaid_epiflu_uploader_v113_surveillance20190109.xls ~/git/CNR-influenza/data/last_gisaid_xls.xls

#files expected at the end of the pipeline
rule all:
    input:
        auspice_jsonH1N1S4 = "auspice/CNR-influenza_H1N1_S4.json",
        auspice_jsonH1N1S6 = "auspice/CNR-influenza_H1N1_S6.json",
        auspice_jsonH3N2S4 = "auspice/CNR-influenza_H3N2_S4.json",
        auspice_jsonH3N2S6 = "auspice/CNR-influenza_H3N2_S6.json",
        auspice_jsonBVICS4 = "auspice/CNR-influenza_BVIC_S4.json",
        auspice_jsonBVICS6 = "auspice/CNR-influenza_BVIC_S6.json",
        convert_csv_ref = "/srv/nfs/ngs-stockage/NGS_Virologie/hcl-vir-ngs/CNRVI/2019_2020/gisaid_epiflu_isolates_H1N1_20200203.xlsx"
        #auspice_jsonBYAMS4 = "auspice/CNR-influenza_BYAM_S4.json",
        #auspice_jsonBYAMS6 = "auspice/CNR-influenza_BYAM_S6.json" 
        #FULL = "auspice/CNR-influenza.json"        



#cp temp_data/ref_data.csv /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/
rule get_last_data_ref:
    input:
        xls_ref_gisaid = "/srv/nfs/ngs-stockage/NGS_Virologie/hcl-vir-ngs/CNRVI/2019_2020/gisaid_epiflu_isolates_H1N1_20200203.xlsx"
    output:
        csv_ref= "temp_data/ref_data.csv"
    shell:
        """
        script/prepare_refxls.py {input} {output} 
        sed -e 's/|/_/g' temp.csv > {output}
        rm temp.csv
       """


rule xls_to_fasta_csv:
    input:
        xls_file = "data/last_gisaid_xls.xls"
    output:
        metadata_raw = "temp_data/metadata_raw.csv",
        fasta_seq = "temp_data/sequences.fasta"      
    shell:
        "script/process_xls.py {input} {output.fasta_seq} {output.metadata_raw}"   

#Produce metadata for all pulled sequences
rule make_metadata_onef:
    input:
        csv_file = rules.xls_to_fasta_csv.output.metadata_raw
    output:
        metadata = "temp_data/metadata.tsv"

    shell:
        "Rscript script/make_meta_onef.R {input} "

#Produce metadata for each subtype and segment
rule make_metadata:
    input:
        csv_file = "temp_data/metadata_raw.csv",
    output:
        H1N1_S4 = "temp_data/H1N1_S4.tsv",
        H1N1_S6 = "temp_data/H1N1_S6.tsv",
        H3N2_S4 = "temp_data/H3N2_S4.tsv",
        H3N2_S6 = "temp_data/H3N2_S6.tsv",
        BVIC_S4 = "temp_data/BVIC_S4.tsv",
        BVIC_S6 = "temp_data/BVIC_S6.tsv",
        BYAM_S4 = "temp_data/BYAM_S4.tsv",
        BYAM_S6 = "temp_data/BYAM_S6.tsv",
               
    shell:
        "Rscript script/make_metadata.R {input.csv_file}  "  

#filter for each dataset using nextstrain-augur
rule augur_filter:
    input:
        seq_file = rules.xls_to_fasta_csv.output.fasta_seq ,
        meta_subtype = "temp_data/{subtype}_{segment}.tsv"

    output:
        filtered_seq = "temp_data/{subtype}_{segment}.fasta"  
    shell:
        "augur filter  "
        "--sequences {input.seq_file} "
        "--metadata {input.meta_subtype}  "
        "--output {output.filtered_seq} " 

#Align each dataset on the reference segment sequence (see readme)
rule augur_align:
    input:
        filter_fasta = rules.augur_filter.output.filtered_seq ,
        ref_seq = "config/{subtype}_{segment}.fasta"
    output:
        align_fasta = "temp_data/{subtype}_{segment}_align.fasta"
    shell:
        "augur align "
        "--sequences {input.filter_fasta} "
        "--reference-sequence {input.ref_seq} "
        "--remove-reference "
        "--fill-gaps "
        "--output {output} "
            
#Build raw tree using augur-nextstrain
rule augur_raw_tree:
    input:
        align_data = rules.augur_align.output.align_fasta
    output:
        raw_tree = "temp_data/{subtype}_{segment}_raw.nwk"
    shell:
        "augur tree "
        "--alignment {input} "
        "--output {output} "

#add temporality to the tree
rule augur_refine:
    input:
        tree = rules.augur_raw_tree.output.raw_tree,
        alignment = rules.augur_align.output.align_fasta,
        meta  = "temp_data/{subtype}_{segment}.tsv"
    output:
        tree = "temp_data/{subtype}_{segment}.nwk",
        node_data = "temp_data/{subtype}_{segment}_branch_lengths.json"
    shell:
        "augur refine "
        "--tree {input.tree} "
        "--alignment {input.alignment} "
        "--metadata {input.meta} "
        "--timetree "
        "--output-tree {output.tree} "
        "--output-node-data {output.node_data} "

#researching nucleotide mutations
rule augur_ancestral:
    input:
        tree = rules.augur_refine.output.tree,
        alignment = rules.augur_align.output.align_fasta
    output:
        node_data = "temp_data/{subtype}_{segment}_nt_muts.json"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data}
        """
#research mutation using genbank data
rule augur_translate:
    input:
        tree = rules.augur_refine.output.tree,
        node_data = rules.augur_ancestral.output.node_data,
        reference = "config/{subtype}_{segment}.gb"
    output:
        node_data = "temp_data/{subtype}_{segment}_aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """
#using clades data from nextstrain on our S4 segments 
rule augur_clades:
    input:
        tree = rules.augur_refine.output.tree,
        aa_muts = rules.augur_translate.output.node_data,
        nuc_muts = rules.augur_ancestral.output.node_data,
        clades = "config/{subtype}_{segment}_clade.tsv"
    output:
        clade_data = "temp_data/{subtype}_{segment}_clades.json"
    wildcard_constraints:
        segment="S4"    
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output {output.clade_data}
        """

#export S4 json for nextstrain visualisation with clades data
rule augur_export_S4:
    input:
        tree = rules.augur_refine.output.tree,
        metadata = "temp_data/{subtype}_{segment}.tsv",
        nt_muts = rules.augur_ancestral.output.node_data,
        aa_muts = rules.augur_translate.output.node_data,
        branch_lengths = rules.augur_refine.output.node_data,
        #clades = rules.augur_clades.output.clade_data
        #auspice_config = "config/auspice_config.json"
    output:
        auspice_json = "auspice/CNR-influenza_{subtype}_{segment}.json",
        #auspice_json = rules.all.input
    wildcard_constraints:
        segment="S4"        
    shell:
        "augur export v2 "
        "--tree {input.tree} "
        "--metadata {input.metadata} "
        "--node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts}  " #{input.clades}
        #"--auspice-config {input.auspice_config} "
        "--output {output.auspice_json} "

#export S6 json for nextstrain visualisation
rule augur_export_S6:
    input:
        tree = rules.augur_refine.output.tree,
        metadata = "temp_data/{subtype}_{segment}.tsv",
        nt_muts = rules.augur_ancestral.output.node_data,
        aa_muts = rules.augur_translate.output.node_data,
        branch_lengths = rules.augur_refine.output.node_data,
        #auspice_config = "config/auspice_config.json"
    output:
        auspice_json = "auspice/CNR-influenza_{subtype}_{segment}.json",
        #auspice_json = rules.all.input
    wildcard_constraints:
        segment="S6"    
    shell:
        "augur export v2 "
        "--tree {input.tree} "
        "--metadata {input.metadata} "
        "--node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts}  " #{input.clades}
        #"--auspice-config {input.auspice_config} "
        "--output {output.auspice_json} "