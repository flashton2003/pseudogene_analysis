from Bio import SeqIO

# Define the input and output file paths
input_file = "GCF_000007545.1.fasta"
output_file = "GCF_000007545.1.frame_shift.fasta"
insertion_position = 399645  # Position to insert 'A'

# Read in the FASTA sequence
with open(input_file, "r") as handle:
    records = list(SeqIO.parse(handle, "fasta"))
    
# Insert 'A' at the specified position in each sequence and write to output
for record in records:
    sequence = record.seq
    modified_sequence = sequence[:insertion_position] + "A" + sequence[insertion_position:]
    record.seq = modified_sequence

# Write the modified sequence to a new FASTA file
with open(output_file, "w") as output_handle:
    SeqIO.write(records, output_handle, "fasta")
