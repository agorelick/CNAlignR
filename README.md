# CNAlign
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

## With the CNalign environment activated, open R, then install additional requirements from bioconductor
```r
# install prerequisites from bioconductor
install.packages("BiocManager")
BiocManager::install("Rsamtools")
BiocManager::install("ggtree")
BiocManager::install("Biobase")

# install ASCAT v3.1.3 from source
install.packages('inst/ascat-3.1.3.tar.gz',type='src',repos=NULL)

# install CNalign R package
install.packages('.',type='src',repos=NULL)
```


# Running CNalign on Whole-Exome Sequencing (WES) data

## Pre-processing input data
For WES input, the input data to CNalign is generated in the same way as for ASCAT: You will use the ASCAT utility prepareHTS() to run _alleleCounter_ (a software which we installed in step XX above) to obtain allele-specific read counts for common SNPs from Tumor/Normal-paired bam files. These will be used as input to CNalign in the following step.

A utility R script is provided in `<CNalign_dir>/scripts/run_ascat_prepareHTS.R`. This can be run from the command line for each pair of tumor/normal bam files for a given patient. For users on HMS's O2 cluster, a template script for the slurm job scheduler is available at `<CNalign_dir>/scripts/run_ascat_prepareHTS_slurm.sh`, which can be used to run run_ascat_prepareHTS.R parallelized with each tumor/normal pair as a distinct submitted job. 

**NB: This Rscript will generate many temporary files into the current directory**, so it's best to `cd` to the intended directory for the output files first, and then run this command as follows. In this example, the directory containing .bam files is one directory up, i.e., the bam file arguments start with ../)

```
Rscript <CNalign_dir>/scripts/run_ascat_prepareHTS.R \
  --patient testpatient \
  --tumor_name Liv3 \
  --tumor_bam ../Liv3_recal.bam \
  --normal_name N1 \
  --normal_bam ../N1_recal.bam \
  --sex XX \
  --genome hg38 \
  --target_bed /n/data1/hms/genetics/naxerova/lab/alex/reference_data/xgen-exome-hyb-panel/xgen-exome-hyb-panel-v2-targets-hg38.bed \
  --allelecounter_exe ~/miniconda3/envs/CNalign/bin/alleleCounter \
  --alleles_prefix /n/data1/hms/genetics/naxerova/lab/alex/reference_data/ascat/G1000_allelesAll_hg38/G1000_alleles_hg38_chr \
  --loci_prefix /n/data1/hms/genetics/naxerova/lab/alex/reference_data/ascat/G1000_lociAll_hg38/G1000_loci_GRCh38_chr
```






