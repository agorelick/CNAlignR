# CNalign
use optimal segment alignment to fit purity and ploidy across multi-region bulk tumor samples

# Installation

## Create a new Conda environment and install dependencies from Conda
conda env create -n "CNalign" -f environment.yml
conda activate CNalign

## in R (inside the CNalign environment) install additional requirements from bioconductor
install.packages("BiocManager")
BiocManager::install("Rsamtools")
BiocManager::install("ggtree")
BiocManager::install("Biobase")

## install ASCAT v3.1.3 from source code in this github repo
install.packges('inst/ascat-3.1.3.tar.gz',type='src',repos=NULL)

