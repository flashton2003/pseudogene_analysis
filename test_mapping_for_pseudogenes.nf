nextflow.enable.dsl=2

// Include original phenix workflow

process PHENIX {
    tag { dataset_id }
    publishDir "${params.root_dir}/${dataset_id}", mode: 'copy'
    container 'flashton/phenix-threaded-samtools:latest'

    input:
    tuple val(dataset_id), path(forward), path(reverse)
    path phenix_config
    path ref

    output:
    path "${dataset_id}*", emit: phenix_output  // Output now within work directory

    script:
    """
    phenix.py prepare_reference -r $ref \
    --mapper bwa \
    --variant gatk
    
    phenix.py run_snp_pipeline \
    -r1 $forward \
    -r2 $reverse \
    -r ${ref} \
    -c ${phenix_config} \
    --keep-temp \
    --sample-name ${dataset_id} \
    -o .  # Save outputs in the current (work) directory

    phenix.py vcf2fasta \
    -i ${dataset_id}.filtered.vcf \
    -o ${dataset_id}_all.fasta \
    --reference ${params.ref} 
    
    gzip ${dataset_id}.vcf
    """
}


process FASTA_TO_FASTQ {
    tag { assembly_id }
    publishDir "${params.root_dir}/${assembly_id}/synthetic_reads", mode: 'copy'
    //container 'quay.io/biocontainers/wgsim:1.0'
    conda 'bioconda::wgsim'

    input:
    tuple val(assembly_id), path(assembly)

    output:
    tuple val(assembly_id), path("${assembly_id}_1.fq"), path("${assembly_id}_2.fq")

    script:
    """
    # Calculate number of reads needed for desired coverage
    GENOME_SIZE=\$(grep -v ">" ${assembly} | tr -d '\\n' | wc -c)
    NUM_READS=\$(echo "scale=0; ${params.coverage} * \$GENOME_SIZE / (2 * ${params.read_length})" | bc)
    
    # Generate paired-end reads using wgsim
    wgsim \
        -1 ${params.read_length} \
        -2 ${params.read_length} \
        -d ${params.fragment_size} \
        -s ${params.fragment_sd} \
        -N \$NUM_READS \
        -e 0 \
        -r 0 \
        ${assembly} \
        ${assembly_id}_1.fq \
        ${assembly_id}_2.fq
    """
}

process ANNOTATE_VCF {
    tag { assembly_id }
    publishDir "${params.root_dir}/${assembly_id}/annotated_vcf", mode: 'copy'
    conda "bioconda::vcf-annotator=0.7"  // Use conda environment for vcf-annotator

    input:
    tuple val(assembly_id), path(vcf)
    path ref_gb  // Reference GenBank file

    output:
    tuple val(assembly_id), path("${assembly_id}.annotated.vcf")

    script:
    """
    vcf-annotator \
        --vcf ${vcf} \
        --genbank ${ref_gb} \
        --output ${assembly_id}.annotated.vcf
    """
}

process IDENTIFY_PSEUDOGENES {
    tag { assembly_id }
    publishDir "${params.root_dir}/${assembly_id}/pseudogenes", mode: 'copy'
    conda "bioconda::pyvcf"

    input:
    tuple val(assembly_id), path(annotated_vcf)

    output:
    tuple val(assembly_id), path("${assembly_id}_pseudogenes.gff"), path("${assembly_id}_pseudogenes.txt")

    script:
    """
    #!/usr/bin/env python3
    import vcf
    import sys
    
    def is_disruptive(variant):
        # Check if variant is likely to disrupt gene function
        effects = variant.INFO.get('ANN', '').split('|')
        disruptive_effects = [
            'frameshift_variant',
            'stop_gained',
            'stop_lost',
            'start_lost',
            'splice_donor_variant',
            'splice_acceptor_variant'
        ]
        return any(effect in effects for effect in disruptive_effects)
    
    # Process VCF file
    vcf_reader = vcf.Reader(filename='${annotated_vcf}')
    pseudogenes = {}
    
    for record in vcf_reader:
        if is_disruptive(record):
            gene_id = record.INFO.get('Gene', 'Unknown')
            if gene_id not in pseudogenes:
                pseudogenes[gene_id] = {
                    'position': record.POS,
                    'variants': []
                }
            pseudogenes[gene_id]['variants'].append(record)
    
    # Write GFF output
    with open('${assembly_id}_pseudogenes.gff', 'w') as gff, \
         open('${assembly_id}_pseudogenes.txt', 'w') as txt:
        
        # Write headers
        gff.write('##gff-version 3\\n')
        txt.write('Gene\\tPosition\\tVariant_Type\\tEffect\\n')
        
        # Write entries
        for gene_id, data in pseudogenes.items():
            for variant in data['variants']:
                # GFF entry
                gff.write(f"{variant.CHROM}\\tpseudogene_pipeline\\tpseudogene\\t" \
                         f"{variant.POS}\\t{variant.POS}\\t.\\t.\\t.\\t" \
                         f"ID={gene_id};variant_type={variant.var_type}\\n")
                
                # Detailed text entry
                txt.write(f"{gene_id}\\t{variant.POS}\\t{variant.var_type}\\t" \
                         f"{variant.INFO.get('ANN', 'Unknown')}\\n")
    """
}

workflow {
    // Create channel from input assemblies
    Channel
        .fromPath(params.assemblies)
        .map { file -> tuple(file.baseName, file) }
        .set { assembly_ch }

    // Generate synthetic FastQ reads
    FASTA_TO_FASTQ(assembly_ch)
        .set { synthetic_reads_ch }

    // Run PHENIX pipeline
    PHENIX(
        synthetic_reads_ch,
        params.phenix_config,
        params.ref
    )
}
