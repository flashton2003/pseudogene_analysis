# Snakefile for genome analysis workflow

# Define the accession and paths directly in the script
GENOME_ACCESSION = ['GCF_000020745.1', 'GCF_000020885.1', 'GCF_000009505.1', 'GCF_000018705.1', 'GCF_000195995.1', 'GCF_000007545.1', 'GCF_000011885.1', 'GCF_000026565.1', 'GCF_000020925.1', 'GCF_000009525.1', 'GCF_000008105.1', 'GCF_000018385.1']
BAKTA_DB = "/data/fast/core/bakta/db"
SALMONELLA_PANGENOME = "/data/fast/salmonella/isangi/pseudogenes/2024.10.03/uniref_taxonomy_id_28901_NOT_name_fr_2024_10_03.fasta"
OUTPUT_DIR = f"/data/fast/salmonella/isangi/pseudogenes/2024.10.29/{GENOME_ACCESSION}"  # Define the output directory variable

rule all:
    input:
        expand(f"{OUTPUT_DIR}/{{accession}}_bakta_db_pseudos.gff", accession=GENOME_ACCESSION),
        expand(f"{OUTPUT_DIR}/{{accession}}_salmonella_pseudos.gff", accession=GENOME_ACCESSION),
        expand(f"{OUTPUT_DIR}/{{accession}}_ncbi_pseudos.gff", accession=GENOME_ACCESSION)

rule download_genome:
    output:
        fasta = f"{OUTPUT_DIR}/raw_data/{{accession}}.fasta",
        protein = f"{OUTPUT_DIR}/raw_data/{{accession}}_protein.faa",
        gff = f"{OUTPUT_DIR}/raw_data/{{accession}}.gff"
    conda:
        "/home/phil/envs/ncbi_datasets.yaml"
    shell:
        """
        datasets download genome accession {wildcards.accession} --include genome,protein,gff3 --filename {wildcards.accession}.zip
        unzip {wildcards.accession}.zip
        mv {wildcards.accession}/data/{wildcards.accession}/*.fna {output.fasta}
        mv {wildcards.accession}/data/{wildcards.accession}/protein.faa {output.protein}
        mv {wildcards.accession}/data/{wildcards.accession}/*.gff {output.gff}
        """

rule annotate_bakta:
    input:
        fasta = f"{OUTPUT_DIR}/raw_data/{{accession}}.fasta"
    output:
        gff = f"{OUTPUT_DIR}/bakta_output/{{accession}}/{{accession}}.gbff"
    conda:
        "/home/phil/envs/bakta.yaml"
    threads: 32
    shell:
        """
        bakta --db {BAKTA_DB} --threads {threads} --output {OUTPUT_DIR}/bakta_output/{wildcards.accession} \
        --prefix {wildcards.accession} --force {input.fasta}
        """

rule run_pseudofinder_bakta_db:
    input:
        gff = rules.annotate_bakta.output.gff
    output:
        gff = f"{OUTPUT_DIR}/{{accession}}_bakta_db_pseudos.gff"
    conda:
        "/home/phil/envs/pseudofinder.yaml"
    threads: 8
    shell:
        """
        python ~/programs/pseudofinder/pseudofinder.py annotate -g {input.gff} \
        --outprefix {OUTPUT_DIR}/{wildcards.accession}_bakta_db \
        --database {BAKTA_DB}/psc.dmnd -di -skpdb -t {threads}
        """

rule run_pseudofinder_salmonella:
    input:
        gff = rules.annotate_bakta.output.gff
    output:
        gff = f"{OUTPUT_DIR}/{{accession}}_salmonella_pseudos.gff"
    conda:
        "/home/phil/envs/pseudofinder.yaml"
    threads: 8
    shell:
        """
        python ~/programs/pseudofinder/pseudofinder.py annotate -g {input.gff} \
        --outprefix {OUTPUT_DIR}/{wildcards.accession}_salmonella \
        --database {SALMONELLA_PANGENOME} -di -t {threads}
        """

rule run_pseudofinder_ncbi:
    input:
        gff = rules.annotate_bakta.output.gff,
        protein = f"{OUTPUT_DIR}/raw_data/{{accession}}_protein.faa"
    output:
        gff = f"{OUTPUT_DIR}/{{accession}}_ncbi_pseudos.gff"
    conda:
        "/home/phil/envs/pseudofinder.yaml"
    threads: 8
    shell:
        """
        python ~/programs/pseudofinder/pseudofinder.py annotate -g {input.gff} \
        --outprefix {OUTPUT_DIR}/{wildcards.accession}_ncbi \
        --database {input.protein} -di -t {threads}
        """

