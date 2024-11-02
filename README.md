# CNalign
use optimal segment alignment to fit purity and ploidy across multi-region bulk tumor samples

# Installation

## Clone the CNalign github repo and cd to it
```
git clone https://github.com/agorelick/CNalign
cd CNalign
```

## Create a new Conda environment and install dependencies from Conda
```
conda env create -n "CNalign" -f environment.yml
conda activate CNalign
```

## In R (inside the CNalign environment) install additional requirements from bioconductor
```r
install.packages("BiocManager")
BiocManager::install("Rsamtools")
BiocManager::install("ggtree")
BiocManager::install("Biobase")
```

## In R (inside the CNalign environment), install ASCAT v3.1.3 from source code in this github repo
```r
install.packges('inst/ascat-3.1.3.tar.gz',type='src',repos=NULL)
```

## In R (inside the CNalign environment), install CNalign R package
```r
devtools::install('.')
```
