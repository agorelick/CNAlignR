#!/bin/bash
#SBATCH -c 20
#SBATCH -t 0-12:00
#SBATCH -p short
#SBATCH --mem 20000
#SBATCH -o run_glimpse_%A_%a.out
#SBATCH --array=1-23
#SBATCH --mail-user=alexander_gorelick@hms.harvard.edu  # Email to which notifications will be sent
#SBATCH --mail-type=ALL

source ~/.bash_profile
mkdir -p glimpse

chr=${SLURM_ARRAY_TASK_ID}
if [[ "$chr" -eq 23 ]]; then
    chr='X'
fi

BAM="original_bams/SD28_NC1_T1.bam"
OUTPUT="glimpse/SD28_NC1"

REF="/home/alg2264/data/alex/reference_data/GLIMPSE_GRCh38/reference_panel/split/1000GP.chr${chr}"
CHUNKS="/home/alg2264/data/alex/reference_data/GLIMPSE_GRCh38/chunks/chunks.chr${chr}.txt"
while IFS="" read -r LINE || [ -n "$LINE" ]; 
do   
    printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)
    CHR=$(echo ${LINE} | cut -d" " -f2)
    REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1)
    REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
    GLIMPSE2_phase_static --bam-file ${BAM} --reference ${REF}_${CHR}_${REGS}_${REGE}.bin --output ${OUTPUT}_${CHR}_${REGS}_${REGE}.bcf --threads 20
done < ${CHUNKS}


