#!/usr/bin/env nextflow

// Enable DSL 2
nextflow.enable.dsl = 2

// Params
params.input_dir = '/data/fast/salmonella/isangi/pseudogenes/2024.11.07'
params.reference_embl = '/data/fast/salmonella/isangi/pseudogenes/2024.11.01/GCF_000022165.1_ASM2216v1_genomic/GCF_000022165.1.embl'
params.reference_faa = '/data/fast/salmonella/isangi/pseudogenes/2024.11.06/2024.11.06.nuccio_baumler_uniprotkb.clean.fasta'
params.nuccio_xlsx = '/data/fast/salmonella/isangi/pseudogenes/2024.11.05/mbo001141769st1.adding_isangi.xlsx'
params.anaerobic_xlsx = '/data/fast/salmonella/isangi/pseudogenes/2024.11.05/mbo001141769st7.central_anaerobic_genes.xlsx'
// params.output_dir = "results"
params.strain_lookup = [
    'GCF_000020705.1': 'SL476',
    'GCF_000020745.1': 'CVM19633',
    'GCF_000020885.1': 'SL483',
    'GCF_000009505.1': 'P125109',
    'GCF_000018705.1': 'SPB7',
    'GCF_000195995.1': 'CT18',
    'GCF_000007545.1': 'Ty2',
    'GCF_000011885.1': 'ATCC 9150',
    'GCF_000020925.1': 'CT_02021853',
    'GCF_000009525.1': '287/91',
    'GCF_000008105.1': 'SC-B67',
    'GCF_000018385.1': 'RKS4594'
]

// Process to run deltaBS.pl in Docker container
process runDeltaBS {
    container 'delta-bit-score'
    
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
        -hd ./ \
        -t /tmp \
        -C 32
    """
}

// Process to run reciprocal diamond analysis
process runReciprocalDiamond {
    conda 'deltaBS_analysis_env.yaml'
    
    input:
    tuple val(strain_id), path(faa_file)
    path(reference_faa)
    
    output:
    tuple val(strain_id), path("${strain_id}_vs_nuccio.reciprocal_diamond.tsv")
    
    script:
    """
    python diamon_reciprocal_best_hits.py \
        --query_fasta ${faa_file} \
        --subject_fasta ${reference_faa} \
        --output ${strain_id}_vs_nuccio.reciprocal_diamond.tsv
    """
}

// Process to join Nuccio and DBS results
process joinNuccioAndDBS {
    conda 'deltaBS_analysis_env.yaml'
    
    input:
    tuple val(strain_id), path(diamond_results), path(dbs_results)
    path(nuccio_xlsx)
    path(anaerobic_xlsx)
    
    output:
    tuple val(strain_id), path("combined_nuccio_and_${strain_id}_dbs.xlsx")
    
    script:
    """
    python join_nuccio_and_dbs.py \
        --nuccio ${nuccio_xlsx} \
        --lookup ${diamond_results} \
        --dbs ${dbs_results} \
        --anaerobic ${anaerobic_xlsx} \
        --output combined_nuccio_and_${strain_id}_dbs.xlsx
    """
}

// Process to calculate summary statistics
process calcSummaryStats {
    conda 'deltaBS_analysis_env.yaml'
    
    input:
    tuple val(strain_id), path(combined_results)
    
    output:
    path("${strain_id}_summary_stats.txt")
    
    script:
    def strain_name = params.strain_lookup[strain_id]
    """
    python calc_deltabs_summary_stats.py \
        ${combined_results} \
        ${strain_name} \
        > ${strain_id}_summary_stats.txt
    """
}

// Main workflow
workflow {
    // Channel for input files
    Channel
        .fromPath("${params.input_dir}/*")
        .filter { it.isDirectory() }
        .map { dir ->
            def strain_id = dir.name
            def embl = file("${dir}/${strain_id}.embl")
            def faa = file("${dir}/${strain_id}.faa")
            return tuple(strain_id, embl, faa)
        }
        .branch {
            embl: it[1].exists()
            faa: it[2].exists()
        }
        .set { input_files }

    // Run deltaBS
    deltaBS_results = runDeltaBS(
        input_files.embl.map { it[0,1] }, 
        params.reference_embl
    )

    // Run reciprocal diamond
    diamond_results = runReciprocalDiamond(
        input_files.faa.map { it[0,2] },
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