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

}
