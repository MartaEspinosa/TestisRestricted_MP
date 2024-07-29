#!/bin/bash

#SBATCH -p long,bigmem            # Partition to submit to
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu 16Gb           # Memory in MB
#SBATCH -J featCounts           # job name
#SBATCH -o /users/genomics/marta/logs/featCounts.%j.out    # File to which standard out will be written
#SBATCH -e /users/genomics/marta/logs/featCounts.%j.err    # File to which standard err will be written

####Change output and output error paths to redirect the log files generated 

###PREPARING NEEDED DATA
PROJECT=$1
DIR=/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/cancers

p=paired-end
strand=unstranded

OUTDIR=$DIR/featureCounts
mkdir -p $OUTDIR

AnnotGTF="/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/human/newReference_Resconstructed/gencode.v38.gffcompare.TestisLiverBrain.annotation.sorted.1transcript.sorted.NOchr.gtf"

module load Subread/2.0.3
########################

# countReadPairs may need to be removed in case of single-end reads
featureCounts -T 10 -p -s 0 -g transcript_id -O --countReadPairs -a $AnnotGTF -o ${OUTDIR}/featureCounts_${PROJECT}.txt /users/genomics/marta/TCGA_RNASeq/${PROJECT}/analysis/05_STAR/uniquely_mapped_2pass_BAM_files/*bam 
# featureCounts -T 10 -p -s 0 -g transcript_id -O --countReadPairs -a $AnnotGTF -o ${OUTDIR}/featureCounts_${PROJECT}.txt /users/genomics/marta/cancers_RNASeq/${PROJECT}/analysis/05_STAR/uniquely_mapped_2pass_BAM_files/*bam 
