library(dplyr)
library(readxl)
library(tidyr)
library(stringr)
library(readr)

# Load Excel file 
dbs_file <- "/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.11.01/typmu_vs_isangi/results.dbs.xlsx" # replace with your actual file path
dbs_data <- read_excel(dbs_file)

# Load GFF file
gff_file <- "/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.11.01/CNS1F3.gff3" # replace with your actual file path
gff_data <- read_delim(gff_file, delim = "\t", comment = "#", col_names = FALSE)

# Rename columns for GFF file
colnames(gff_data) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# Extract `locus_tag` from `attributes` column in GFF file
gff_data <- gff_data %>%
  filter(feature == "CDS") %>%
  mutate(locus_tag = str_extract(attributes, "locus_tag=[^;]+")) %>%
  mutate(locus_tag = str_remove(locus_tag, "locus_tag="))

# Join Excel data with GFF data based on gene_2 and locus_tag
merged_data <- dbs_data %>%
  left_join(gff_data, by = c("gene_2" = "locus_tag"))

# Split the `attributes` column into multiple columns by `;`
merged_data <- merged_data %>%
  separate(attributes, into = paste0("attribute_", 1:10), sep = ";", fill = "right")

# Save the merged data to a new Excel file
write.csv(merged_data, "/Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.11.01/typmu_vs_isangi/results.dbs.annot.csv", row.names = FALSE)
