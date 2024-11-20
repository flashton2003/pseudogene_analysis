# Pseudogene Detection Pipeline

This repository contains a collection of scripts and workflows for detecting and analyzing pseudogenes in bacterial genomes, specifically focusing on Salmonella strains. It integrates multiple tools and approaches to perform comprehensive pseudogene analysis, allowing comparison and validation of results from different detection methods.

### Features

The pipeline supports three primary pseudogene detection approaches:

Bakta Annotation: Genome annotation using Bakta, with pseudogene detection.
PseudoFinder: Pseudogene detection tailored to specific bacterial genomes.
Delta-Bitscore (DBS) Analysis: Uses DIAMOND reciprocal best hits to calculate delta-bitscores for identifying pseudogenes.
It also provides utilities for combining and validating results from these methods, generating statistics, and producing visualizations.

### Overall workflow

1a.pseudofinder_bakta_workflow.smk runs pseudofinder and bakta
1b.deltaBS_workflow.nf runs deltaBS

Then, there are two ways those results are assessed against truth i) by matching the coordinates of each call to the coordinates of the calls in the truth set, and ii) by matching the protein sequence of each CDS against the truth set.

2a.genomic_coords_join_validation_with_nuccio.py matches based on co-ordinates.
2b.0.diamond_best_hits.py does the diamond best hit analysis.
2b.1.diamond_join_validation_with_nuccio.py - joins the results to the truth using the diamond best hits.

Then, stats

3.pseudogene_stats.py

Plotting

4.pseudogene_stat_plotting.py

and finally, if you want to run the diamond workflow on a non-truth set genome

5.diamond_join_with_nuccio.py

Excel files: Combined results, summary statistics, and annotations.
TSV files: Intermediate results and processed datasets.
Plots (PNG): Visualizations of comparative statistics, such as sensitivity vs. PPV.

### Command line examples

`for sample in GCF_000007545.1 GCF_000008105.1 GCF_000009505.1 GCF_000009525.1 GCF_000011885.1 GCF_000018385.1 GCF_000018705.1 GCF_000020705.1 GCF_000020745.1 GCF_000020885.1 GCF_000020925.1 GCF_000026565.1 GCF_000195995.1; do python scripts/2b.0.diamond_best_hits.py --query_fasta ./2024.10.29/$sample/bakta_output/$sample.faa --subject_fasta 2024.11.06/2024.11.06.nuccio_baumler_uniprotkb.clean.fasta --output 2024.11.14/${sample}_vs_nuccio.diamond.tsv; done`

`for sample in GCF_000007545.1 GCF_000008105.1 GCF_000009505.1 GCF_000009525.1 GCF_000011885.1 GCF_000018385.1 GCF_000018705.1 GCF_000020705.1 GCF_000020745.1 GCF_000020885.1 GCF_000020925.1 GCF_000026565.1 GCF_000195995.1; do python scripts/2b.1.diamond_join_validation_with_nuccio.py  --nuccio 2024.11.05b/mbo001141769st1.adding_isangi.xlsx --bakta 2024.10.29/pseudogene_calls/${sample}.bakta.gff3 --pseudofinder-baktadb 2024.10.29/pseudogene_calls/${sample}_bakta_db_pseudos.gff --pseudofinder-salmonella 2024.10.29/pseudogene_calls/${sample}_salmonella_pseudos.gff --pseudofinder-ncbi 2024.10.29/pseudogene_calls/${sample}_ncbi_pseudos.gff  --diamond 2024.11.14/${sample}_vs_nuccio.diamond.tsv --anaerobic 2024.11.05b/mbo001141769st7.central_anaerobic_genes.xlsx --dbs 2024.11.07/$sample/results.dbs --output 2024.11.14/$sample.calls_vs_nuccio.xlsx; done`

`for sample in GCF_000007545.1 GCF_000008105.1 GCF_000009505.1 GCF_000009525.1 GCF_000011885.1 GCF_000018385.1 GCF_000018705.1 GCF_000020705.1 GCF_000020745.1 GCF_000020885.1 GCF_000020925.1 GCF_000026565.1 GCF_000195995.1; do python scripts/2b.1.diamond_join_validation_with_nuccio.py  --nuccio 2024.11.05b/mbo001141769st1.adding_isangi.xlsx --bakta 2024.10.29/pseudogene_calls/$sample.bakta.gff3 --pseudofinder-baktadb 2024.10.29/pseudogene_calls/${sample}_bakta_db_pseudos.gff --pseudofinder-salmonella 2024.10.29/pseudogene_calls/${sample}_salmonella_pseudos.gff --pseudofinder-ncbi 2024.10.29/pseudogene_calls/${sample}_ncbi_pseudos.gff  --anaerobic 2024.11.05b/mbo001141769st7.central_anaerobic_genes.xlsx --dbs 2024.11.07/$sample/results.dbs --output 2024.11.14b/$sample.calls_vs_nuccio.coords.xlsx --gcf $sample; done`

`python scripts/3.pseudogene_stats.py --input_dir 2024.11.14/ --output_file 2024.11.14/2024.11.14.pseudogene_validation_results.inc_other_groups.diamond.csv`

`python scripts/4.pseudogene_stat_plotting.py --input_file 2024.11.14b/2024.11.14.pseudogene_validation_results.inc_other_groups.coords.csv --ppv_plot 2024.11.14b/2024.11.14.pseudogene_validation_results.coords.sens_vs_ppv.png --truth_plot 2024.11.14b/2024.11.14.pseudogene_validation_results.coords.truth_vs_total_calls.png --cam_plot 2024.11.14b/2024.11.14.pseudogene_validation_results.coords.cam_truth_vs_calls.png`

`python scripts/F1.diamond_join_with_nuccio.py --nuccio 2024.11.05b/mbo001141769st1.adding_isangi.xlsx --pseudofinder-baktadb 2024.11.14c/CIV13RE2_pseudofinder_pseudos.gff --diamond 2024.11.14c/CIV13RE2_vs_nuccio.diamond.tsv --output 2024.11.14c/CIV13RE2_pf_baktadb_vs_nuccio.xlsx --anaerobic 2024.11.05b/mbo001141769st7.central_anaerobic_genes.xlsx`

### Input files

[mbo001141769st7.central_anaerobic_genes.xlsx
]([url](https://www.dropbox.com/scl/fi/r5ftw4kb99g96bds5q0hx/mbo001141769st7.central_anaerobic_genes.xlsx?rlkey=6pg5z9vd1gnye2gzo6x5fnkfn&dl=0))[mbo001141769st1.adding_isangi.xlsx]([url](https://www.dropbox.com/scl/fi/d4u0x5lnki8j7huq9nhzm/mbo001141769st1.adding_isangi.xlsx?rlkey=gsqermq1dlvpgpawy3uauyttf&dl=0))
[2024.11.06.nuccio_baumler_uniprotkb.clean.fasta]([url](https://www.dropbox.com/scl/fi/2kayz9u194zhut9d57xli/2024.11.06.nuccio_baumler_uniprotkb.clean.fasta?rlkey=65f9y17i13gcu7um30dwdhenm&dl=0))
