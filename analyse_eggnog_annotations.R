# Load required libraries
library(dplyr)
library(clusterProfiler)
library(KEGGREST)  # For accessing KEGG pathways directly

# Function to read and process eggNOG mapper output
read_eggnog <- function(file_path) {
  # Read eggNOG mapper output, tab-delimited
  eggnog_data <- read.delim(file_path, 
                            comment.char = "#", 
                            stringsAsFactors = FALSE)
  
  # Rename columns to standard eggNOG mapper output
  colnames(eggnog_data)[1:3] <- c("query", "seed_ortholog", "evalue")
  
  # Extract KEGG annotations from the KEGG column
  # Assuming KEGG annotations are in column 12 (adjust if needed)
  kegg_data <- eggnog_data %>%
    select(query, KEGG = 12) %>%
    filter(!is.na(KEGG) & KEGG != "") %>%
    # Split multiple KEGG annotations
    tidyr::separate_rows(KEGG, sep = ",")
  
  # Clean KEGG IDs (remove 'ko:' prefix if present)
  kegg_data$KEGG <- gsub("ko:", "", kegg_data$KEGG)
  
  return(kegg_data)
}

# Function to get KEGG pathway information for prokaryotes
get_kegg_pathways <- function() {
  # Get all KEGG pathways
  pathways <- keggList("pathway")
  
  # Get pathway to KO mapping
  path2ko <- lapply(names(pathways), function(pid) {
    tryCatch({
      path <- keggGet(pid)[[1]]
      if (!is.null(path$GENE)) {
        kos <- path$GENE[grep("^K[0-9]", names(path$GENE))]
        data.frame(
          pathway_id = pid,
          pathway_name = pathways[pid],
          ko = names(kos),
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) NULL)
  })
  
  # Combine all pathways
  path2ko_df <- do.call(rbind, path2ko)
  
  # Clean pathway IDs
  path2ko_df$pathway_id <- gsub("path:", "", path2ko_df$pathway_id)
  path2ko_df$ko <- gsub("ko:", "", path2ko_df$ko)
  
  return(path2ko_df)
}

# Function to perform KEGG enrichment analysis
perform_kegg_enrichment <- function(gene_subset, background_genes, kegg_data, pathway_info,
                                    p_value = 0.05, q_value = 0.2) {
  # Create gene to KEGG mapping
  gene2kegg <- kegg_data %>%
    filter(query %in% background_genes) %>%
    select(KEGG, query) %>%
    distinct()
  
  # Prepare input for enrichment analysis
  gene_list <- gene_subset
  universe <- background_genes
  
  # Perform enrichment analysis
  enrich_result <- enricher(gene = gene_list,
                            universe = universe,
                            TERM2GENE = gene2kegg,
                            TERM2NAME = pathway_info %>%
                              select(pathway_id, pathway_name) %>%
                              distinct(),
                            pvalueCutoff = p_value,
                            qvalueCutoff = q_value)
  
  # Convert results to data frame
  result_df <- as.data.frame(enrich_result)
  
  return(result_df)
}


# Example usage:
# First get the pathway information (this step takes time)
# You only need to run this once and can save the results
message("Retrieving KEGG pathways... This may take several minutes.")
pathway_info <- get_kegg_pathways()
# Optionally save for future use
saveRDS(pathway_info, "kegg_pathway_info.rds")
# For future runs, you can load with:
# pathway_info <- readRDS("kegg_pathway_info.rds")

# Replace 'path/to/eggnog_output.tsv' with your actual file path
file_path <- "path/to/eggnog_output.tsv"

# Read and process eggNOG mapper output
kegg_annotations <- read_eggnog(file_path)

# Example: create a subset of genes (replace with your actual gene subset)
# This could be your genes of interest
all_genes <- unique(kegg_annotations$query)
gene_subset <- sample(all_genes, size = round(length(all_genes) * 0.1))  # 10% random sample for example

# Perform enrichment analysis using the retrieved pathway info
enrichment_results <- perform_kegg_enrichment(
  gene_subset = gene_subset,
  background_genes = all_genes,
  kegg_data = kegg_annotations,
  pathway_info = pathway_info
)

# View results
print(enrichment_results)

# Create visualization
if (nrow(enrichment_results) > 0) {
  # Create dot plot of enriched pathways
  library(ggplot2)
  
  ggplot(enrichment_results, 
         aes(x = Count, y = reorder(Description, Count))) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8)) +  # Adjust text size for readability
    labs(x = "Gene Count",
         y = "KEGG Pathway",
         title = "KEGG Pathway Enrichment Analysis in Salmonella",
         subtitle = paste("Total genes analyzed:", length(gene_subset)),
         color = "Adjusted\np-value",
         size = "Gene\nCount")
  
  # Save plot if needed
  # ggsave("salmonella_kegg_enrichment_plot.pdf", width = 12, height = 8)
}

# Export results to CSV
# write.csv(enrichment_results, "salmonella_kegg_enrichment_results.csv", row.names = FALSE)