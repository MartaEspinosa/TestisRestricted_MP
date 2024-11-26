#!/bin/bash

#SBATCH -p long,bigmem            # Partition to submit to
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu 15Gb           # Memory in MB
#SBATCH -J bam2fastq_Wang           # job name
#SBATCH -o /users/genomics/marta/logs/bam2fastq_Wang_%A_%a.out    # File to which standard out will be written
#SBATCH -e /users/genomics/marta/logs/bam2fastq_Wang_%A_%a.err    # File to which standard err will be written

INDIR=/users/genomics/saraa/projectTestis/fastq/RnaSeq

declare -a sample_array

for file in $INDIR/human*fastq.gz; do
    sample=${file%%_r1*}
    sample=${sample##*/}
    
    sample_array=(${sample_array[@]} $sample) # Add new element at the end of the array
done

i=$(($SLURM_ARRAY_TASK_ID -1))

sample=${sample_array[i]}

FQDIR=/users/genomics/saraa/projectTestis/fastq/RnaSeq

#### ---------------- STAR --------------- ###
module load STAR/2.7.8a-GCC-10.2.0

BAMDIR=/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/human/STAR

######################################################################################################
#####################################ALIGNMENT########################################################
GENOMEDIR=/data/genomics/marta/genomes

ANNOTGENE=$GENOMEDIR/Annot_files_GTF/GENCODE_47_M36
GNMIDX=$GENOMEDIR/Index_Genomes_STAR_2.8/Idx_Gencode_v38_hg38_readlength75

#gzip fastq files are considered in the code, as well as paired-end reads samples.
#two pass mode is activated
#output with uniquely mapped reads ONLY
STAR --runThreadN $SLURM_CPUS_PER_TASK --limitBAMsortRAM 50000000000 \
    --genomeDir $GNMIDX --readFilesCommand zcat --readFilesIn $FQDIR/${sample}_r1.fastq.gz --outFileNamePrefix\
    ${BAMDIR}/$sample --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outFilterType BySJout\
    --outFilterMultimapNmax 1 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999\
    --outFilterMismatchNoverLmax 0.05 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
    --sjdbGTFfile $ANNOTGENE/gencode.v47.primary_assembly.annotation.fixed.gtf --twopassMode Basic

######################################################################################################
#########################################INDEX########################################################

module purge
module load SAMtools/1.12-GCC-10.2.0

samtools index ${BAMDIR}/${sample}Aligned.sortedByCoord.out.bam ${BAMDIR}/${sample}Aligned.sortedByCoord.out.bai