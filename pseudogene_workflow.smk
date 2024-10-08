# Snakefile for genome analysis workflow

# Define the accession and paths directly in the script
GENOME_ACCESSION = "GCF_000006945.2"
BAKTA_DB = "/data/fast/core/bakta/db"
SALMONELLA_PANGENOME = "/data/fast/salmonella/isangi/pseudogenes/2024.10.03/uniref_taxonomy_id_28901_NOT_name_fr_2024_10_03.fasta"

rule all:
    input:
        expand("{accession}_bakta_db_pseudos.gff", accession=GENOME_ACCESSION),
        expand("{accession}_salmonella_pseudos.gff", accession=GENOME_ACCESSION),
        expand("{accession}_ncbi_pseudos.gff", accession=GENOME_ACCESSION)

rule download_genome:
    output:
        fasta = "raw_data/{accession}.fasta",
        protein = "raw_data/{accession}_protein.faa",
        gff = "raw_data/{accession}.gff"
    conda:
        "/home/phil/envs/ncbi_datasets.yaml"
    shell:
        """
        datasets download genome accession {wildcards.accession} --include genome,protein,gff3
        unzip ncbi_dataset.zip
        mv ncbi_dataset/data/{wildcards.accession}/*.fna {output.fasta}
        mv ncbi_dataset/data/{wildcards.accession}/protein.faa {output.protein}
        mv ncbi_dataset/data/{wildcards.accession}/*.gff {output.gff}
        rm -rf ncbi_dataset ncbi_dataset.zip
        """

rule annotate_bakta:
    input:
        fasta = "raw_data/{accession}.fasta"
    output:
        gff = "bakta_output/{accession}/{accession}.gff3"
    conda:
        "/home/phil/envs/bakta.yaml"
    threads: 32
    shell:
        """
        bakta --db {BAKTA_DB} --threads {threads} --output bakta_output/{wildcards.accession} \
        --prefix {wildcards.accession} --force {input.fasta}
        """

rule run_pseudofinder_bakta_db:
    input:
        gbff = "bakta_output/{accession}/{accession}.gff"
    output:
        gff = "{accession}_bakta_db_pseudos.gff"
    conda:
        "/home/phil/envs/pseudofinder.yaml"
    threads: 8
    shell:
        """
        python ~/programs/pseudofinder/pseudofinder.py annotate -g {input.gff} \
        --outprefix {wildcards.accession}_bakta_db \
        --database {BAKTA_DB}/psc.dmnd -di -skpdb -t {threads}
        mv {wildcards.accession}_bakta_db_pseudos.gff {output.gff}
        """

rule run_pseudofinder_salmonella:
    input:
        gff = "bakta_output/{accession}/{accession}.gff"
    output:
        gff = "{accession}_salmonella_pseudos.gff"
    conda:
        "/home/phil/envs/pseudofinder.yaml"
    threads: 8
    shell:
        """
        python ~/programs/pseudofinder/pseudofinder.py annotate -g {input.gff} \
        --outprefix {wildcards.accession}_salmonella \
        --database {SALMONELLA_PANGENOME} -di -t {threads}
        mv {wildcards.accession}_salmonella_pseudos.gff {output.gff}
        """

rule run_pseudofinder_ncbi:
    input:
        gff = "bakta_output/{accession}/{accession}.gff",
        protein = "raw_data/{accession}_protein.faa"
    output:
        gff = "{accession}_ncbi_pseudos.gff"
    conda:
        "/home/phil/envs/pseudofinder.yaml"
    threads: 8
    shell:
        """
        python ~/programs/pseudofinder/pseudofinder.py annotate -g {input.gff} \
        --outprefix {wildcards.accession}_ncbi \
        --database {input.protein} -di -t {threads}
        mv {wildcards.accession}_ncbi_pseudos.gff {output.gff}
        """
