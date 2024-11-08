#library(readr)
library(tidyverse)
library(clusterProfiler)
#library(enrichplot)
library(data.table)
library(KEGGREST)

# protein ids of HDCs from here - /Users/flashton/Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.11.01/typmu_vs_isangi/results.dbs.xlsx
protein_ids <- c('BIANFB_12320', 'BIANFB_14145', 'BIANFB_19365', 'BIANFB_14920', 'BIANFB_14920', 'BIANFB_17960', 'BIANFB_22230', 'BIANFB_17415', 'BIANFB_20545', 'BIANFB_18930', 'BIANFB_18485', 'BIANFB_02390', 'BIANFB_21655', 'BIANFB_21500', 'BIANFB_20515', 'BIANFB_20515', 'BIANFB_11945', 'BIANFB_13485', 'BIANFB_19280', 'BIANFB_06665', 'BIANFB_01830', 'BIANFB_16760', 'BIANFB_19725', 'BIANFB_15575', 'BIANFB_13220', 'BIANFB_21490', 'BIANFB_20990', 'BIANFB_20440', 'BIANFB_19100', 'BIANFB_17895', 'BIANFB_13730', 'BIANFB_08715', 'BIANFB_16900', 'BIANFB_21960', 'BIANFB_21960', 'BIANFB_14680', 'BIANFB_09340', 'BIANFB_03955', 'BIANFB_19715', 'BIANFB_21340', 'BIANFB_18155', 'BIANFB_19985', 'BIANFB_19020', 'BIANFB_22280', 'BIANFB_10685', 'BIANFB_10685', 'BIANFB_19975', 'BIANFB_06890', 'BIANFB_21930', 'BIANFB_16800', 'BIANFB_22935', 'BIANFB_19540', 'BIANFB_08920', 'BIANFB_08920', 'BIANFB_00390', 'BIANFB_07105', 'BIANFB_13550', 'BIANFB_18520', 'BIANFB_19295', 'BIANFB_21460', 'BIANFB_21190', 'BIANFB_20460', 'BIANFB_14070', 'BIANFB_18085', 'BIANFB_07645', 'BIANFB_20890', 'BIANFB_12355', 'BIANFB_22185', 'BIANFB_03255', 'BIANFB_04190', 'BIANFB_15610', 'BIANFB_01830', 'BIANFB_16530', 'BIANFB_01830', 'BIANFB_03295', 'BIANFB_08400', 'BIANFB_06665', 'BIANFB_21275', 'BIANFB_21775', 'BIANFB_22135', 'BIANFB_20810', 'BIANFB_00345', 'BIANFB_17715', 'BIANFB_21245', 'BIANFB_08880', 'BIANFB_21630', 'BIANFB_16415', 'BIANFB_14885', 'BIANFB_03785', 'BIANFB_17345', 'BIANFB_16640', 'BIANFB_03285', 'BIANFB_03175', 'BIANFB_16185', 'BIANFB_01960', 'BIANFB_22805', 'BIANFB_14150', 'BIANFB_14150', 'BIANFB_22560', 'BIANFB_07655', 'BIANFB_13660', 'BIANFB_08135', 'BIANFB_07635', 'BIANFB_06665', 'BIANFB_03250', 'BIANFB_02455', 'BIANFB_16830', 'BIANFB_09500', 'BIANFB_08515', 'BIANFB_02010', 'BIANFB_10260', 'BIANFB_02815', 'BIANFB_03275', 'BIANFB_03985', 'BIANFB_06340', 'BIANFB_02960', 'BIANFB_05410', 'BIANFB_10315', 'BIANFB_04810', 'BIANFB_09730', 'BIANFB_12445', 'BIANFB_10495', 'BIANFB_19665', 'BIANFB_02010', 'BIANFB_23275', 'BIANFB_21740', 'BIANFB_06850', 'BIANFB_15190', 'BIANFB_16790', 'BIANFB_12390', 'BIANFB_09330', 'BIANFB_01960', 'BIANFB_03450', 'BIANFB_07680', 'BIANFB_09115', 'BIANFB_11590', 'BIANFB_15920', 'BIANFB_01585', 'BIANFB_08560', 'BIANFB_17765', 'BIANFB_02685', 'BIANFB_07245', 'BIANFB_03990', 'BIANFB_02605', 'BIANFB_20930', 'BIANFB_08460', 'BIANFB_04640', 'BIANFB_07655', 'BIANFB_06275', 'BIANFB_18635', 'BIANFB_00750', 'BIANFB_15735', 'BIANFB_19875', 'BIANFB_01175', 'BIANFB_16125', 'BIANFB_08665', 'BIANFB_10635', 'BIANFB_15125', 'BIANFB_10010', 'BIANFB_03250', 'BIANFB_11115', 'BIANFB_16400', 'BIANFB_08860', 'BIANFB_17390', 'BIANFB_11265', 'BIANFB_09625', 'BIANFB_18625', 'BIANFB_00875', 'BIANFB_00655', 'BIANFB_04920', 'BIANFB_20825', 'BIANFB_20945', 'BIANFB_19035', 'BIANFB_19035', 'BIANFB_13370', 'BIANFB_09700', 'BIANFB_05195', 'BIANFB_21965', 'BIANFB_08735', 'BIANFB_09005', 'BIANFB_15760', 'BIANFB_08905', 'BIANFB_00850', 'BIANFB_17770', 'BIANFB_14280', 'BIANFB_10830', 'BIANFB_20155', 'BIANFB_22330', 'BIANFB_18335', 'BIANFB_18335', 'BIANFB_18335', 'BIANFB_08855', 'BIANFB_13285', 'BIANFB_06090', 'BIANFB_14230', 'BIANFB_01935', 'BIANFB_05380', 'BIANFB_15460', 'BIANFB_07110', 'BIANFB_11290', 'BIANFB_06565', 'BIANFB_16690', 'BIANFB_00395', 'BIANFB_08020', 'BIANFB_20070', 'BIANFB_09875', 'BIANFB_05140')

# read the data

CNS1F3_eggnog_emapper_annotations <- read_delim("Dropbox/GordonGroup/STRATAA_XDR_Salmonella_Isangi/pseudogene_finding/2024.11.05/CNS1F3_eggnog.emapper.annotations.tsv", 
                                                delim = "\t", escape_double = FALSE, 
                                                comment = "##", trim_ws = TRUE)
CNS1F3_eggnog_emapper_annotations <- CNS1F3_eggnog_emapper_annotations %>% rename(query = `#query`)

## this section uses the clusterProfiler package

#get columns 1 (query) and 9 (KO terms)
kegg_data <- CNS1F3_eggnog_emapper_annotations[c(1,13)]

# expand, since some genes/proteins will have multiple assigned KO terms
kegg_data <- data.table(kegg_data)
kegg_data <- kegg_data[, list(KEGG_Pathway = unlist(strsplit(KEGG_Pathway , ","))), by = query]

print(length(unique(kegg_data$KEGG_Pathway)))

# select the needed columns
kegg_data <- kegg_data[,c(2,1)]

#kegg_data_filtered <- kegg_data %>% filter(KEGG_ko != '-')
#protein_ids_filterd <- protein_ids[protein_ids %in% kegg_data_filtered$query]

enr_res <- enricher(protein_ids, TERM2GENE=kegg_data, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, minGSSize = 10)

#enr_res <- enricher(protein_ids, TERM2GENE=kegg_data, pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 10)


write.table(enr_res, file = '~/Desktop/tmp.csv')


## i also want to do a quick validation, because i'm not quite sure what's going on with cluster profiler.

# Calculate pathway counts for genes of interest and overall
enrichment_results <- kegg_data %>%
  mutate(is_interest = ifelse(query %in% protein_ids, "interest", "other")) %>%
  group_by(KEGG_Pathway, is_interest) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = is_interest, values_from = count, values_fill = 0) %>%
  rename(in_interest = interest, in_other = other)

# Run chi-square test for each pathway
enrichment_results <- enrichment_results %>%
  rowwise() %>%
  mutate(
    chi_sq_pval = chisq.test(matrix(c(in_interest, in_other,
                                      sum(in_interest), sum(in_other)), 
                                    nrow = 2))$p.value
  ) %>%
  ungroup()

# Adjust p-values for multiple testing using Benjamini-Hochberg (FDR)
enrichment_results <- enrichment_results %>%
  mutate(FDR = p.adjust(chi_sq_pval, method = "BH"))

# Select and arrange the output
enrichment_results <- enrichment_results %>%
  select(KEGG_pathway, in_interest, in_other, chi_sq_pval, FDR) %>%
  arrange(FDR)

## some exploratory code to look up more info about the kegg pathway from the ID using the rest API.

# Function to get parent terms for a KO ID
get_ko_hierarchy <- function(ko_id) {
  # Add 'ko:' prefix if not present
  if(!grepl("^ko:", ko_id)) {
    ko_id <- paste0("ko:", ko_id)
  }
  
  # Get KO info
  ko_info <- keggGet(ko_id)[[1]]
  
  # Extract BRITE hierarchy information
  brite_info <- ko_info$BRITE
  
  # Process and clean up the hierarchy information
  hierarchy_levels <- lapply(brite_info, function(x) {
    # Split by newlines and remove empty strings
    lines <- strsplit(x, "\n")[[1]]
    lines <- lines[lines != ""]
    
    # Remove the KO line itself (last line)
    lines <- lines[-length(lines)]
    
    # Clean up the formatting
    lines <- gsub("^[A-Z] +", "", lines)
    
    return(lines)
  })
  
  return(hierarchy_levels)
}

# Get hierarchy for K04023
hierarchy <- get_ko_hierarchy("ko:K04023")
ko_info <- keggGet("ko:K04023")
# Print results
print(ko_info[[1]]$BRITE[length(ko_info[[1]]$BRITE) - 1])
