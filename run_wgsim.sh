# Define variables
assembly="GCF_000007545.1.frame_shift.fasta"    # Path to the genome assembly file
assembly_id="GCF_000007545.1.frame_shift"            # Identifier for the output files
params_coverage=30                      # Desired coverage
params_read_length=150                  # Read length for each read
params_fragment_size=500                # Mean fragment size
params_fragment_sd=50                   # Fragment size standard deviation

# Calculate genome size
GENOME_SIZE=$(grep -v ">" ${assembly} | tr -d '\n' | wc -c)

# Calculate the number of reads needed for the specified coverage
NUM_READS=$(echo "scale=0; ${params_coverage} * $GENOME_SIZE / (2 * ${params_read_length})" | bc)

# Generate paired-end reads using wgsim
wgsim \
    -1 ${params_read_length} \
    -2 ${params_read_length} \
    -d ${params_fragment_size} \
    -s ${params_fragment_sd} \
    -N $NUM_READS \
    -e 0 \
    -r 0 \
    ${assembly} \
    ${assembly_id}_1.fq \
    ${assembly_id}_2.fq

