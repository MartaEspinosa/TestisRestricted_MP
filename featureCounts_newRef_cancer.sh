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
DIR=/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/cancers

p=$2
strand=$3

OUTDIR=$DIR/featureCounts
mkdir -p $OUTDIR

# AnnotGTF="/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/human/newReference_Resconstructed/gencode.v38.gffcompare.TestisLiverBrain.annotation.sorted.1transcript.sorted.NOchr.gtf"

AnnotGTF="/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/human/newReference_Resconstructed/gencode.v47.gffcompare.TestisLiverBrain.annotation.sorted.1transcript.sorted.NOchr.fixed.gtf"
module load Subread/2.0.3
########################

# DIR=/users/genomics/marta/TCGA_RNASeq
# DIR=/users/genomics/marta/cancers_RNASeq
# DIR=/projects_eg/projects/marta
DIR=/users/genomics/marta

if [ $p == "paired-end" ]; then
    if [ $strand == "firststrand" ]; then
        featureCounts -T 10 -p -s 2 -g transcript_id -O --countReadPairs -a $AnnotGTF -o ${OUTDIR}/featureCounts_${PROJECT}.txt ${DIR}/${PROJECT}/analysis/05_STAR/uniquely_mapped_2pass_BAM_files/*bam 
    elif [ $strand == "secondstrand" ]; then
        featureCounts -T 10 -p -s 1 -g transcript_id -O --countReadPairs -a $AnnotGTF -o ${OUTDIR}/featureCounts_${PROJECT}.txt ${DIR}/${PROJECT}/analysis/05_STAR/uniquely_mapped_2pass_BAM_files/*bam 
    elif [ $strand == "unstranded" ]; then
        featureCounts -T 10 -p -s 0 -g transcript_id -O --countReadPairs -a $AnnotGTF -o ${OUTDIR}/featureCounts_${PROJECT}.txt ${DIR}/${PROJECT}/analysis/05_STAR/uniquely_mapped_2pass_BAM_files/*bam 
    fi
fi
if [ $p == "single-end" ]; then
    if [ $strand == "firststrand" ]; then
        featureCounts -T 10 -s 2 -g transcript_id -O -a $AnnotGTF -o ${OUTDIR}/featureCounts_${PROJECT}.txt ${DIR}/${PROJECT}/analysis/05_STAR/uniquely_mapped_2pass_BAM_files/*bam 
    elif [ $strand == "secondstrand" ]; then
        featureCounts -T 10 -s 1 -g transcript_id -O -a $AnnotGTF -o ${OUTDIR}/featureCounts_${PROJECT}.txt ${DIR}/${PROJECT}/analysis/05_STAR/uniquely_mapped_2pass_BAM_files/*bam 
    elif [ $strand == "unstranded" ]; then
        featureCounts -T 10 -s 0 -g transcript_id -O -a $AnnotGTF -o ${OUTDIR}/featureCounts_${PROJECT}.txt ${DIR}/${PROJECT}/analysis/05_STAR/uniquely_mapped_2pass_BAM_files/*bam 
    fi
fi