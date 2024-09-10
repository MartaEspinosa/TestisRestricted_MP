#!/bin/bash

#SBATCH -p long,bigmem            # Partition to submit to
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu 15Gb           # Memory in MB
#SBATCH -J bam2fastq_Wang           # job name
#SBATCH -o /users/genomics/marta/logs/bam2fastq_Wang_%A_%a.out    # File to which standard out will be written
#SBATCH -e /users/genomics/marta/logs/bam2fastq_Wang_%A_%a.err    # File to which standard err will be written

specie=$1
# INDIR=/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/$specie/RiboSeq/bam
INDIR=/datasets/marta/Wang2020/RiboSeq
declare -a sample_array

for file in $INDIR/${specie}*bam; do
    sample=${file%%.bam*}
    sample=${sample##*/}

    sample_array=(${sample_array[@]} $sample) # Add new element at the end of the array
done

i=$(($SLURM_ARRAY_TASK_ID -1))

sample=${sample_array[i]}

# INDIR=/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/$specie/RiboSeq/bam
INDIR=/datasets/marta/Wang2020/RiboSeq
FQDIR=/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/$specie/RiboSeq/fastq
mkdir -p $FQDIR

#### ---------------- BAM 2 FASTQ --------------- ###
module load Java/11.0.2
module load picard/2.25.1-Java-11

java -jar $EBROOTPICARD/picard.jar SamToFastq I=$INDIR/${sample}.bam FASTQ=$FQDIR/${sample}_r1.fastq
# module load biobambam2
#
# mkdir -o $FQDIR/unmatched
#
# bamtofastq filename=$bam F=$FQDIR/${sample}_r1.fastq O=$FQDIR/unmatched/${sample}_unmatched_r1.fastq

