#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-1:30
#SBATCH -p short
#SBATCH --mem 32000
#SBATCH --array=0-10 
#SBATCH -o slurm/run_ascat_prepareHTS_%A_%a.out
#SBATCH --mail-user=alexander_gorelick@hms.harvard.edu  # Email to which notifications will be sent
#SBATCH --mail-type=ALL

# create directory for slurm log files
mkdir -p slurm 

# patient-specific info
declare -a tumor_samples=("PM2" "PM3a" "PM3b" "PM4a" "PM4b" "PT1P1" "PT2P3" "PT3P3" "PT5P2" "PT6P5" "PT7P1") # tumor samples only
patient="LM31"      # patient ID
normal_sample="N2"  # normal sample
sex="XX"            # sex (XX or XY)
bam_dir="/n/data1/hms/genetics/naxerova/lab/alex/peritoneal_revision_WES/FC_08873/preprocessing_ds" # directory with .bam files

# get the tumor bam/name from the SLURM array
i=${SLURM_ARRAY_TASK_ID}
tumor_sample=${tumor_samples[$i]}

# get full paths for given pair of T/N bam files
tumor_bam=${bam_dir}"/${tumor_sample}_recal.bam"
normal_bam=${bam_dir}"/${normal_sample}_recal.bam"

# run ASCAT/alleleCounter to get input data for CNalign
Rscript <CNalign_dir>/scripts/run_ascat_prepareHTS.R \
    --patient $patient \
    --tumor_name $tumor_sample \
    --tumor_bam $tumor_bam \
    --normal_name $normal_sample \
    --normal_bam $normal_bam \
    --sex $sex \
    --genome $genome \
    --target_bed /n/data1/hms/genetics/naxerova/lab/alex/reference_data/xgen-exome-hyb-panel/xgen-exome-hyb-panel-v2-targets-hg38.bed \
    --allelecounter_exe ~/miniconda3/envs/CNalign/bin/alleleCounter \
    --alleles_prefix /n/data1/hms/genetics/naxerova/lab/alex/reference_data/ascat/G1000_allelesAll_hg38/G1000_alleles_hg38_chr \
    --loci_prefix /n/data1/hms/genetics/naxerova/lab/alex/reference_data/ascat/G1000_lociAll_hg38/G1000_loci_GRCh38_chr



