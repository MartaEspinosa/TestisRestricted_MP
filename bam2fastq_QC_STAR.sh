#!/bin/bash

#SBATCH -p long,bigmem            # Partition to submit to
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu 15Gb           # Memory in MB
#SBATCH -J bam2fastq_Wang           # job name
#SBATCH -o /users/genomics/marta/logs/bam2fastq_Wang_%A_%a.out    # File to which standard out will be written
#SBATCH -e /users/genomics/marta/logs/bam2fastq_Wang_%A_%a.err    # File to which standard err will be written

INDIR=/projects_eg/projects/marta/Wang2020_download_testisCellTypes/RNASeq

declare -a sample_array

for file in $INDIR/*bam; do
    sample=${file%%.bam*}
    sample=${sample##*/}
    
    sample_array=(${sample_array[@]} $sample) # Add new element at the end of the array
done

i=$(($SLURM_ARRAY_TASK_ID -1))

sample=${sample_array[i]}

INDIR=/projects_eg/projects/marta/Wang2020_download_testisCellTypes/RNASeq
FQDIR=/users/genomics/marta/Wang2020_download_testisCellTypes/RNASeq/fastq/
mkdir -p $FQDIR

QCDIR=/users/genomics/marta/Wang2020_download_testisCellTypes/RNASeq/QC
mkdir -p $QCDIR

#### ---------------- BAM 2 FASTQ --------------- ###
module load Java/11.0.2
module load picard/2.25.1-Java-11

java -jar $EBROOTPICARD/picard.jar SamToFastq I=$INDIR/${sample}.bam FASTQ=$FQDIR/${sample}_r1.fastq


#### ---------------- QC --------------- ###
module purge  ## Why? Clear out .bashrc /.bash_profile settings that might interfere
module load FastQ-Screen/0.14.1
module load Bowtie2/2.4.2-GCC-10.2.0             # Required for Fastqscreen

mkdir -p $QCDIR/fastqscreen 
config=/datasets/FastQ_Screen_Genomes/fastq_screen_cova.conf #para mouse y macaque

echo "Starting FastQScreen"
fastq_screen --threads $SLURM_CPUS_PER_TASK --conf $config --outdir $QCDIR/fastqscreen $FQDIR/${sample}_r1.fastq
echo "FastQScreen done"


module load FastQC/0.11.7-Java-1.8.0_162
mkdir -p $QCDIR/fastqc

echo "Starting FastQC"
fastqc -t $SLURM_CPUS_PER_TASK -I -O $QCDIR/fastqc $FQDIR/${sample}_r1.fastq 
echo "FastQC done"

#### ---------------- STAR --------------- ###
module load STAR/2.7.8a-GCC-10.2.0

BAMDIR=/users/genomics/marta/Wang2020_download_testisCellTypes/RNASeq/STAR

######################################################################################################
#####################################ALIGNMENT########################################################
GENOMEDIR=/users/genomics/saraa/projectTestis/genomes/mouse

ANNOTGENE=$GENOMEDIR
GNMIDX=$GENOMEDIR/Index_Genomes_STAR_2.8/Idx_GRCm39_readlength75

#gzip fastq files are considered in the code, as well as paired-end reads samples.
#two pass mode is activated
#output with uniquely mapped reads ONLY
STAR --runThreadN $SLURM_CPUS_PER_TASK --limitBAMsortRAM 50000000000 \
    --genomeDir $GNMIDX --readFilesIn $FQDIR/${sample}_r1.fastq --outFileNamePrefix\
    ${BAMDIR}/$sample --outSAMattributes All --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outFilterType BySJout\
    --outFilterMultimapNmax 1 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999\
    --outFilterMismatchNoverLmax 0.05 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000\
    --sjdbGTFfile /datasets/Common/ReferenceGenomes/Mouse/mm39/GenomesAndAnnotation/gencode.vM32.annotation.gtf --twopassMode Basic

######################################################################################################
#########################################INDEX########################################################

module purge
module load SAMtools/1.12-GCC-10.2.0

samtools index ${BAMDIR}/${sample}Aligned.sortedByCoord.out.bam ${BAMDIR}/${sample}Aligned.sortedByCoord.out.bai