#!/bin/bash

dir=/datasets/marta/ATACseq/Corces2018/bigWigs

projects="BRCA BLCA LIHC COAD PRAD LUSC LUAD KIRC"

module load ucsc_tools

for proj in $projects; do
  echo $proj

  ## bigWigs to bedGraph
  bigWigMerge $dir/$proj/*bw $dir/$proj/${proj}.bedGraph

  ## Sort bedGraph
  sort -k1,1 -2,2n $dir/$proj/${proj}.bedGraph > $dir/$proj/${proj}.sorted.bedGraph

  ## bedGraphToBigWig
  bedGraphToBigWig $dir/$proj/${proj}.sorted.bedGraph https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes $dir/$proj/${proj}.bigWig

done
