#!/bin/bash
#SBATCH -c 4
#SBATCH -t 0-0:20
#SBATCH -p short
#SBATCH --mem 32000
#SBATCH -o get_CNalign_obj.out
#SBATCH --mail-user=alexander_gorelick@hms.harvard.edu  # Email to which notifications will be sent
#SBATCH --mail-type=ALL

# load R with CNalign installed
module load R/4.3.1
Rscript get_CNalign_obj.R

