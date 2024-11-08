#!/usr/bin/env nextflow

// Enable DSL 2
nextflow.enable.dsl = 2

// Process to run deltaBS.pl in Docker container
process runDeltaBS {
    container 'delta-bit-score'
    containerOptions = '-v /data/fast/salmonella/isangi/pseudogenes/2024.11.01:/mnt/deltaBS'
    
    input:
    tuple val(strain_id), path(embl_file)
    path(reference_embl)
    
    output:
    tuple val(strain_id), path("${strain_id}/results.dbs")
    
    script:
    """
    mkdir -p ${strain_id}
    ~/deltaBS.pl -f embl \
        -f1 ${reference_embl} \
        -f2 ${embl_file} \
        -o ${strain_id} \
        -hp /usr/bin/ \
        -hd /mnt/deltaBS \
        -t /tmp \
        -C 32
   """
}

// Process to run reciprocal diamond analysis
process runReciprocalDiamond {
    conda 'scripts/deltaBS_analysis_env.yml'
    
    input:
    tuple val(strain_id), path(faa_file)
    path(reference_faa)
    
    output:
    tuple val(strain_id), path("${strain_id}_vs_nuccio.reciprocal_diamond.tsv")
    
    script:
    """
    python /data/fast/salmonella/isangi/pseudogenes/scripts/diamon_reciprocal_best_hits.py \
        --query_fasta ${faa_file} \
        --subject_fasta ${reference_faa} \
        --output ${strain_id}_vs_nuccio.reciprocal_diamond.tsv
    """
}

// Process to join Nuccio and DBS results
process joinNuccioAndDBS {
    conda 'scripts/deltaBS_analysis_env.yml'
    
    input:
    tuple val(strain_id), path(diamond_results), path(dbs_results)
    path(nuccio_xlsx)
    path(anaerobic_xlsx)
    
    output:
    tuple val(strain_id), path("combined_nuccio_and_${strain_id}_dbs.xlsx")
    
    script:
    """
    python scripts/join_nuccio_and_dbs.py \
        --nuccio ${nuccio_xlsx} \
        --lookup ${diamond_results} \
        --dbs ${dbs_results} \
        --anaerobic ${anaerobic_xlsx} \
        --output combined_nuccio_and_${strain_id}_dbs.xlsx
    """
}

// Process to calculate summary statistics
process calcSummaryStats {
    conda 'scripts/deltaBS_analysis_env.yml'
    
    input:
    tuple val(strain_id), path(combined_results)
    
    output:
    path("${strain_id}_summary_stats.txt")
    
    script:
    def strain_name = params.strain_lookup[strain_id]
    """
    python scripts/calc_deltabs_summary_stats.py \
        ${combined_results} \
        ${strain_name} \
        > ${strain_id}_summary_stats.txt
    """
}

// Main workflow
workflow {
    // Create channels for EMBL and FAA files
    embl_files = Channel
        .fromPath("${params.input_dir}/*.embl")
        .map { file -> 
            def strain_id = file.name.toString().tokenize('.')[0]
            return tuple(strain_id, file)
        }

    faa_files = Channel
        .fromPath("${params.input_dir}/*.faa")
        .map { file -> 
            def strain_id = file.name.toString().tokenize('.')[0]
            return tuple(strain_id, file)
        }

    // Run deltaBS
    deltaBS_results = runDeltaBS(
        embl_files,
        params.reference_embl
    )

    // Run reciprocal diamond
    diamond_results = runReciprocalDiamond(
        faa_files,
        params.reference_faa
    )

    // Join the results
    combined_results = joinNuccioAndDBS(
        deltaBS_results.join(diamond_results),
        params.nuccio_xlsx,
        params.anaerobic_xlsx
    )

    // Calculate summary statistics
    calcSummaryStats(combined_results)
}
