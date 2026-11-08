# Assignment-4-Introduction-to-R
# DNA Sequence Analysis in R (chr1_GL383518v1_alt)

## Identifying Information
* **Programmer:** Qaliya Omar
* **Language:** R (Base R, Version 4.0 or higher)
* **Date Submitted:** November 08, 2026
* **Purpose:** This script reads a biological FASTA sequence file, performs string indexing, generates a Watson-Crick-Franklin reverse complement sequence, calculates nucleotide frequency distributions across 1000 base pair windows using lists, and constructs data frames to evaluate sequence block sums.

## Required Files
* **`chr1_GL383518v1_alt.fa`**: FASTA formatted text file containing the sequence for Human Chromosome 1 secondary assembly `GL383518v1`. Must be stored in the same working directory as `dna_analysis.R`.

## Software & Library Requirements
* **Language Environment:** Standard R installation
* **Required Libraries:** None. Built entirely using base R functions (no external packages required).

## Instructions for Execution
1. Open your R console, RStudio, or terminal command prompt.
2. Set your working directory to the folder containing `dna_analysis.R` and `chr1_GL383518v1_alt.fa`:
   ```R
   setwd("/path/to/your/folder")
3. Run the R script:
    ```R
   source("dna_analysis.R")
