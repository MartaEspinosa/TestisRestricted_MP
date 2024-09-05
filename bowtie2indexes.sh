#!/bin/bash

#SBATCH -p long,bigmem            # Partition to submit to
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu 15Gb           # Memory in MB
#SBATCH -J bowtie2_index           # job name
#SBATCH -o /users/genomics/marta/logs/bowtie2_index_%j.out    # File to which standard out will be written
#SBATCH -e /users/genomics/marta/logs/bowtie2_index_%j.err    # File to which standard err will be written

module load Bowtie2

bowtie2-build -f /datasets/Common/ReferenceGenomes/Chicken/GRCg7b/GenomesAndAnnotation/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa /users/genomics/marta/bowtie2_indexes/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b
