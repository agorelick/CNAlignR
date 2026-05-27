#!/bin/bash
#SBATCH -c 2
#SBATCH -t 0-8:00
#SBATCH -p short
#SBATCH --mem 8000
#SBATCH -o run_QDNAseq.out
#SBATCH --mail-user=alexander_gorelick@hms.harvard.edu  # Email to which notifications will be sent
#SBATCH --mail-type=ALL

module load R/4.3.1
Rscript run_QDNAseq.R
