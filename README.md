# CNAlign
use optimal segment alignment to fit purity and ploidy across multi-region bulk tumor samples


## 1. Installation

### Clone the CNAlignR github repo and cd to it
```bash
git clone https://github.com/agorelick/CNAlignR
cd CNAlignR
```

### Create a new Conda environment and install dependencies from Conda
```bash
conda env create -n "CNAlignR" -f environment.yml
conda activate CNAlignR
```

### With the CNAlignR environment activated, open R, then install additional requirements from bioconductor
```r
# install prerequisites from bioconductor
install.packages("BiocManager")
BiocManager::install("Rsamtools")
BiocManager::install("ggtree")
BiocManager::install("Biobase")

# install ASCAT v3.1.3 from source
install.packages('inst/ascat-3.1.3.tar.gz',type='src',repos=NULL)

# install CNAlignR R package
install.packages('.',type='src',repos=NULL)
```


## 2. Running CNAlignR on Whole-Exome Sequencing (WES) data

### 2.1 Pre-processing input data
For WES input, the input data to CNAlignR is generated in the same way as for ASCAT: You will use the ASCAT utility prepareHTS() to run _alleleCounter_ (a software which we installed in step XX above) to obtain allele-specific read counts for common SNPs from Tumor/Normal-paired bam files. These will be used as input to CNAlignR in the following step.

A utility R script is provided in `<CNAlignR_dir>/scripts/run_ascat_prepareHTS.R`. This can be run from the command line for each pair of tumor/normal bam files for a given patient. For users on HMS's O2 cluster, a template script for the slurm job scheduler is available at `<CNAlignR_dir>/scripts/run_ascat_prepareHTS_slurm.sh`, which can be used to run run_ascat_prepareHTS.R parallelized with each tumor/normal pair as a distinct submitted job. 

**NB: This Rscript will generate many temporary files into the current directory**, so it's best to `cd` to the intended directory for the output files first, and then run this command as follows. In this example, the directory containing .bam files is one directory up, i.e., the bam file arguments start with ../). 

**NB:** In the following command, encode patient sex as "XX" or "XY" (this will be used for the expected copies in chrX/chrY). Also, encode the reference genome version/build as "hg19" or "hg38", regardless of whether the .bam files have "chr"-prefixed chromosome names.

```bash
Rscript <CNAlignR_dir>/scripts/run_ascat_prepareHTS.R \
  --patient testpatient \
  --tumor_name Liv3 \
  --tumor_bam ../Liv3_recal.bam \
  --normal_name N1 \
  --normal_bam ../N1_recal.bam \
  --sex XX \
  --genome hg38 \
  --target_bed /n/data1/hms/genetics/naxerova/lab/alex/reference_data/xgen-exome-hyb-panel/xgen-exome-hyb-panel-v2-targets-hg38.bed \
  --allelecounter_exe ~/miniconda3/envs/CNAlignR/bin/alleleCounter \
  --alleles_prefix /n/data1/hms/genetics/naxerova/lab/alex/reference_data/ascat/G1000_allelesAll_hg38/G1000_alleles_hg38_chr \
  --loci_prefix /n/data1/hms/genetics/naxerova/lab/alex/reference_data/ascat/G1000_lociAll_hg38/G1000_loci_GRCh38_chr
```

Once this has finished running, your directory should be populated with **six files for each tumor sample**, with names in the format:
```bash
<PatientID>_<TumorID>_<NormalID>_Tumor_LogR.txt
<PatientID>_<TumorID>_<NormalID>_Tumor_BAF.txt
<PatientID>_<TumorID>_<NormalID>_Tumor_BAF_rawBAF.txt
<PatientID>_<TumorID>_<NormalID>_Germline_LogR.txt
<PatientID>_<TumorID>_<NormalID>_Germline_BAF.txt
<PatientID>_<TumorID>_<NormalID>_Germline_BAF_rawBAF.txt
```

### 2.2 Generate CNAlignR input data (.rds files)

With the CNAlignR conda environment still activated, use the provided R script to merge these sample files together and create the input data R object for CNAlignR. This will likely require **at least 32GB of RAM** and may take **a few hours to complete** (depending on the number of samples).
```bash
# example command to generate input data object for one patient.
Rscript ~/repos/CNAlignR/scripts/merge_alleleCounter_data.R \
  --patient C66         # <PatientID> \
  --normal_sample c66N3 # <NormalID>  \
  --sex XX              # XX|XY       \
  --build hg19          # hg19|hg38   \
  --GCcontentfile ~/data/alex/reference_data/ascat/GC_G1000_hg19.txt \
  --replictimingfile ~/data/alex/reference_data/ascat/RT_G1000_hg19.txt \
  --obj_file            # <PatientID>_CNAlignR_obj.rds
```

Your directory should now have three .rds files with the (default) filenames as below. **Consider these the input data for CNAlignR.**
```bash
# CNAlignR input with high-penalty mpcf (multi-sample segment alignment, penalty=300).
# This should have fewer copy number segments (ideally < 150), making it easier to fit purity/ploidy in each sample.
<PatientID>_CNAlignR_obj_mpcf.rds         

# CNAlignR input with low-penalty mpcf (penalty=25).
# This may have hundreds of copy number segments, making it more sensitive for focal copy number changes. Use these segments for assigning copy number to specific genes/mutations, based on the purity obtained from the aforementioned high-penalty version.
<PatientID>_CNAlignR_obj_mpcf_hisens.rds

# same as the above outputs but without running MPCF. This is mainly for QC/recovery in case MPCF fails due to excessive time/memory usage. 
<PatientID>_CNAlignR_obj.rds              
```



