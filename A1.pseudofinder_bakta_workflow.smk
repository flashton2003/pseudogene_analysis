import os

# Define the accession and paths directly in the script

# Get single accession from config
GENOME_ACCESSION = config.get("accession", "")
if not GENOME_ACCESSION:
    raise ValueError("Please provide a genome accession using --config accession=GCF_XXXXXX")


BAKTA_DB = "/data/fast/core/bakta/db"
SALMONELLA_PANGENOME = "/data/fast/salmonella/isangi/pseudogenes/2024.10.03/uniref_taxonomy_id_28901_NOT_name_fr_2024_10_03.fasta"
OUTPUT_BASE = "/data/fast/salmonella/isangi/pseudogenes/2024.10.29"

rule all:
    input:
        expand(f"{OUTPUT_BASE}/{{accession}}/{{accession}}_bakta_db_pseudos.gff", accession=GENOME_ACCESSION),
        expand(f"{OUTPUT_BASE}/{{accession}}/{{accession}}_salmonella_pseudos.gff", accession=GENOME_ACCESSION),
        expand(f"{OUTPUT_BASE}/{{accession}}/{{accession}}_ncbi_pseudos.gff", accession=GENOME_ACCESSION)

rule download_genome:
    output:
        fasta = f"{OUTPUT_BASE}/{{accession}}/raw_data/{{accession}}.fasta",
        protein = f"{OUTPUT_BASE}/{{accession}}/raw_data/{{accession}}_protein.faa",
        gff = f"{OUTPUT_BASE}/{{accession}}/raw_data/{{accession}}.gff"
    conda:
        "/home/phil/envs/ncbi_datasets.yaml"
    shell:
        """
        datasets download genome accession {wildcards.accession} --include genome,protein,gff3 --filename {wildcards.accession}.zip
        unzip {wildcards.accession}.zip -d {wildcards.accession}
        mv {wildcards.accession}/ncbi_dataset/data/{wildcards.accession}/*.fna {output.fasta}
        mv {wildcards.accession}/ncbi_dataset/data/{wildcards.accession}/protein.faa {output.protein}
        mv {wildcards.accession}/ncbi_dataset/data/{wildcards.accession}/*.gff {output.gff}
        """

rule annotate_bakta:
    input:
        fasta = f"{OUTPUT_BASE}/{{accession}}/raw_data/{{accession}}.fasta"
    output:
        gff = f"{OUTPUT_BASE}/{{accession}}/bakta_output/{{accession}}.gbff"
    conda:
        "/home/phil/envs/bakta.yaml"
    threads: 32
    shell:
        """
        bakta --db {BAKTA_DB} --threads {threads} --output {OUTPUT_BASE}/{wildcards.accession}/bakta_output \
        --prefix {wildcards.accession} --force {input.fasta}
        """

rule run_pseudofinder_bakta_db:
    input:
        gff = rules.annotate_bakta.output.gff
    output:
        gff = f"{OUTPUT_BASE}/{{accession}}/{{accession}}_bakta_db_pseudos.gff"
    conda:
        "/home/phil/envs/pseudofinder.yaml"
    threads: 8
    shell:
        """
        python ~/programs/pseudofinder/pseudofinder.py annotate -g {input.gff} \
        --outprefix {OUTPUT_BASE}/{wildcards.accession}/{wildcards.accession}_bakta_db \
        --database {BAKTA_DB}/psc.dmnd -di -skpdb -t {threads}
        """

rule run_pseudofinder_salmonella:
    input:
        gff = rules.annotate_bakta.output.gff
    output:
        gff = f"{OUTPUT_BASE}/{{accession}}/{{accession}}_salmonella_pseudos.gff"
    conda:
        "/home/phil/envs/pseudofinder.yaml"
    threads: 8
    shell:
        """
        python ~/programs/pseudofinder/pseudofinder.py annotate -g {input.gff} \
        --outprefix {OUTPUT_BASE}/{wildcards.accession}/{wildcards.accession}_salmonella \
        --database {SALMONELLA_PANGENOME} -di -t {threads}
        """

rule run_pseudofinder_ncbi:
    input:
        gff = rules.annotate_bakta.output.gff,
        protein = f"{OUTPUT_BASE}/{{accession}}/raw_data/{{accession}}_protein.faa"
    output:
        gff = f"{OUTPUT_BASE}/{{accession}}/{{accession}}_ncbi_pseudos.gff"
    conda:
        "/home/phil/envs/pseudofinder.yaml"
    threads: 8
    shell:
        """
        python ~/programs/pseudofinder/pseudofinder.py annotate -g {input.gff} \
        --outprefix {OUTPUT_BASE}/{wildcards.accession}/{wildcards.accession}_ncbi \
        --database {input.protein} -di -t {threads}
        """

