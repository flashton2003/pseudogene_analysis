from Bio import SeqIO
from collections import defaultdict

# Input and output file names
input_file = "/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.11.06/2024.11.06.nuccio_baumler_uniprotkb.fasta"
output_file = "/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.11.06/2024.11.06.nuccio_baumler_uniprotkb.clean.fasta"

# Dictionary to store unique sequences with their ID
sequences = defaultdict(list)

# Process the input FASTA file
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        # Check if the ID contains '|' and should be simplified
        if "|" in record.id:
            id_part = record.id.split('|')[1]  # Extract the middle part of the ID
        else:
            id_part = record.id  # Keep the original ID if no '|' is present

        # Check if the sequence already exists for this ID
        if any(rec.seq == record.seq for rec in sequences[id_part]):
            continue  # Skip duplicate sequence
        elif id_part in sequences and any(rec.seq != record.seq for rec in sequences[id_part]):
            print(f"Warning: ID {id_part} has differing sequences.")
        
        # Add the record with the processed ID to the dictionary
        record.id = id_part
        record.description = ""
        sequences[id_part].append(record)

    # Write unique records to the output file
    SeqIO.write((rec for rec_list in sequences.values() for rec in rec_list), outfile, "fasta")