# Pseudogene Detection Pipeline

This is an LLM generated README, it's generally sensible, but interpret with some caution.

This repository contains a collection of scripts for detecting and analyzing pseudogenes in bacterial genomes, specifically focusing on Salmonella strains. The pipeline implements three different approaches to pseudogene detection:

1. PseudoFinder and Bakta-based detection
2. Delta-bitscore (DBS) analysis
3. Combined analysis of both approaches

## Requirements

### Environment Setup
A Conda environment file (`deltaBS_analysis_env.yml`) is provided with the required dependencies:

```yaml
dependencies:
  - diamond
  - python>=3.8
  - pandas>=1.5.0
  - numpy>=1.21.0
  - seaborn>=0.12.0
  - matplotlib>=3.6.0
  - openpyxl
```

To create the environment:
```bash
conda env create -f deltaBS_analysis_env.yml
```

## Scripts Overview

### A. PseudoFinder and Bakta Analysis
0. `A1.pseudofinder_bakta_workflow.smk`: runs pseudofinder and bakta on the nuccio and baumler genome set.
1. `A2.join_nuccio_and_pseudofinder_or_bakta.py`: Combines Nuccio dataset with PseudoFinder or Bakta results
2. `A3.pseudofinder_bakta_summary_stats.py`: Generates summary statistics for PseudoFinder and Bakta predictions

### B. Delta-bitscore (DBS) Analysis
0. `B1.deltaBS_workflow.nf`: runs deltaBS.
1. `B2.diamon_reciprocal_best_hits.py`: Performs DIAMOND reciprocal best hits analysis to enable matching between my annotated genome and the nuccio & baumler "[proteome]([url](https://www.dropbox.com/scl/fi/da92gzzesypw0tax46qsm/2024.11.06.nuccio_baumler_uniprotkb.clean.fasta?rlkey=zk2wtryx17q3mh0ugndf1fjxj&dl=0))".
2. `B3.join_nuccio_and_dbs.py`: Combines Nuccio dataset with DBS results
3. `B4.deltaBS_summary_stats.py`: Generates summary statistics for DBS predictions
4. `B5.deltaBS_overall_plot.py`: Creates visualization plots for DBS analysis

### C. Combined Analysis
1. `C1.combine_bakta_and_non_bakta.py`: Integrates results from both Bakta and non-Bakta approaches
2. `C2.plotting_pseudogene_results.py`: Generates comparative plots of pseudogene detection results
3. `C3.plotting_pseudogene_combined_bakta_results.py`: Creates visualizations for combined analysis results

## Usage

### 1. PseudoFinder and Bakta Analysis
```bash
# Join Nuccio and PseudoFinder/Bakta results
python A2.join_nuccio_and_pseudofinder_or_bakta.py --truth <truth_file> \
    --calls <calls_file> --output <output_file>

# Generate summary statistics
python A3.pseudofinder_bakta_summary_stats.py <excel_file> <reference_name>
```

### 2. Delta-bitscore Analysis
```bash
# Run DIAMOND reciprocal best hits
python B2.diamon_reciprocal_best_hits.py --query_fasta <query.faa> \
    --subject_fasta <subject.faa> --output <output.tsv>

# Join Nuccio and DBS results
python B3.join_nuccio_and_dbs.py --nuccio <nuccio_file> \
    --lookup <lookup_file> --dbs <dbs_file> \
    --anaerobic <anaerobic_file> --output <output_file>

# Generate DBS statistics
python B4.deltaBS_summary_stats.py <excel_file> <reference_name>
```

### 3. Combined Analysis
```bash
# Combine Bakta and non-Bakta results
python C1.combine_bakta_and_non_bakta.py <bakta_file> <nonbakta_file> <output_dir>

# Generate plots
python C2.plotting_pseudogene_results.py
python C3.plotting_pseudogene_combined_bakta_results.py
```

## Configuration

The repository includes two configuration files:
- `deltaBS_nextflow.config`: Configuration for DBS analysis pipeline
- `nextflow.config`: General Nextflow configuration

## Output

The scripts generate various output files including:
- Excel spreadsheets with combined analysis results
- TSV files with summary statistics
- PNG files with visualization plots
- Comparative analysis results between different methods

## Notes

- The pipeline is specifically designed for Salmonella strains and includes mappings for various reference genomes
- The analysis takes into account both extraintestinal (EI) and gastrointestinal (GI) Salmonella strains
- Results include sensitivity and positive predictive value (PPV) calculations for each method

