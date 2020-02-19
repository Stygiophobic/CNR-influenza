#singularity shell /srv/nfs/ngs-stockage/NGS_Virologie/NEXTSTRAIN/nextstrainV4.simg
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
        #auspice_jsonBYAMS4 = "auspice/CNR-influenza_BYAM_S4.json",
        #auspice_jsonBYAMS6 = "auspice/CNR-influenza_BYAM_S6.json" 
        #FULL = "auspice/CNR-influenza.json" 
    params:
        data_rep = "/srv/nfs/ngs-stockage/NGS_Virologie/hcl-vir-ngs/CNRVI/2019_2020/"               



#Get gisaid data made by Gwendo
rule get_last_data_gisaid:
    output:
        xls_file = "data/last_gisaid_xls.xls" 
    shell:
        """
        cp {rules.all.params.data_rep}gisaid_epiflu_uploader*[0-9].xls data/
        mv data/gisaid_epiflu_uploader*[0-9].xls data/last_gisaid_xls.xls
        """        

rule xls_to_fasta_csv:
    input:
        xls_file = rules.get_last_data_gisaid.output.xls_file
    output:
        metadata_raw = "temp_data/metadata_raw.csv",
        fasta_seq = "temp_data/sequences.fasta"      
    shell:
        "script/process_xls.py {input} {output.fasta_seq} {output.metadata_raw}"   

rule process_ref_data:
    output:
        ref_H1N1 = "temp_data/meta_H1N1.csv",
        ref_H3N2 = "temp_data/meta_H3N2.csv",
        ref_BVIC = "temp_data/meta_BVIC.csv",
        ref_BYAM = "temp_data/meta_BYAM.csv",
        ref_seq = "temp_data/SEQ_REF.fasta"
    shell:
        """
        sed -e 's/|/_/g' {rules.all.params.data_rep}REF_DATA_NEXTSTRAIN/SEQREF_H1N1.csv > {output.ref_H1N1}
        sed -e 's/|/_/g' {rules.all.params.data_rep}REF_DATA_NEXTSTRAIN/SEQREF_H3N2.csv > {output.ref_H3N2}
        sed -e 's/|/_/g' {rules.all.params.data_rep}REF_DATA_NEXTSTRAIN/SEQREF_BVIC.csv > {output.ref_BVIC}
        sed -e 's/|/_/g' {rules.all.params.data_rep}REF_DATA_NEXTSTRAIN/SEQREF_BYAM.csv > {output.ref_BYAM}
        cat {rules.all.params.data_rep}REF_DATA_NEXTSTRAIN/SEQREF_*.fasta > {output.ref_seq}
        """        


#Produce metadata for each subtype and segment
rule make_metadata:
    input:
        csv_file = "temp_data/metadata_raw.csv",
        ref_H1N1 = rules.process_ref_data.output.ref_H1N1,
        ref_H3N2 = rules.process_ref_data.output.ref_H3N2,
        ref_BVIC = rules.process_ref_data.output.ref_BVIC,
        ref_BYAM = rules.process_ref_data.output.ref_BYAM ,   

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
        "Rscript script/make_metadata.R {input.csv_file} {input.ref_H1N1} "
        "{input.ref_H3N2} {input.ref_BVIC} {input.ref_BYAM} "  

rule merge_fasta:
    input:
        ref_seq = rules.process_ref_data.output.ref_seq,
        data_seq = rules.xls_to_fasta_csv.output.fasta_seq
    output:
        all_seq = "temp_data/ALL_SEQ_REF.fasta"   
    shell:
        """
        cat {input.ref_seq} {input.data_seq}  > temp_data/temp_merge.fasta 
        sed -e 's/|/_/g' temp_data/temp_merge.fasta > temp_data/temp_merge2.fasta
        tool/seqkit rmdup --by-name --ignore-case temp_data/temp_merge2.fasta > {output}  
        """    

#filter for each dataset using nextstrain-augur
rule augur_filter:
    input:
        seq_file = rules.merge_fasta.output.all_seq ,
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