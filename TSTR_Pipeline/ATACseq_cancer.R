# Version D With ChIPpeakAnno
library("ChIPpeakAnno")
library("ChIPseeker")
library("tidyr")
library("dplyr")
library("GenomicFeatures")
library("BSgenome.Hsapiens.UCSC.hg38")

# Get chromosome lengths
chr_lengths <- as.data.frame(seqlengths(Hsapiens))
chr_lengths_mb <- chr_lengths / 1e6  # Convert base pairs to megabases

chr_lengths_mb <- chr_lengths_mb %>% subset(!grepl("_", rownames(chr_lengths_mb)))
chr_lengths_mb$seqnames <- rownames(chr_lengths_mb)
names(chr_lengths_mb)[1] <- "length"

# Load annotation
gtf_file <- "/data/genomics/marta/genomes/Annot_files_GTF/GENCODE_47_M36/gencode.v47.primary_assembly.annotation.fixed.gtf"
txdb     <- makeTxDbFromGFF(gtf_file)
annotatio_granges <- transcripts(txdb)
names(annotatio_granges) <- mcols(annotatio_granges)$tx_name

# Read ATAC files
atac.peak.files <- c(
  list.files(path="/datasets/marta/ATACseq/Corces2018", all.files=T, full.names = T,pattern = "bed") ## PANCANCER PEAKS 
)

peaks_annotation <- NULL
for (peak_file in atac.peak.files){
  print(peak_file)
  peaks <- readPeakFile(peak_file)
  peaks_df <- as.data.frame(peaks)
    
  ## peaks by length
  peaksXchr <- as.data.frame(table(peaks_df$seqnames))
  names(peaksXchr) <- c("seqnames","num_peaks")
  peaksXchr <- merge(chr_lengths_mb, peaksXchr, by="seqnames")

  peaksXchr$peakDensity <- peaksXchr$num_peaks / peaksXchr$length
  
  ann.tmp <- annotatePeakInBatch(myPeakList = peaks,AnnotationData = annotatio_granges ,output = "upstream&inside",multiple = T)
  peaks_annotation  <- c(peaks_annotation , ann.tmp)
}

peaks_annotation <- do.call(c, peaks_annotation)
peaks_annotation_df <- as.data.frame(peaks_annotation)
peaks_annotation_df$feature <- gsub("\\..","",peaks_annotation_df$feature)


### Candidates 
tumorReactDIR = "/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/cancers/log2ratio3x/cancertypes"
candidatesORFs = read.csv(file.path(tumorReactDIR,"TSTR_candidatesORFs_fullcharacterized.csv"))
genes_candidatesORFs = candidatesORFs[,c("transcript_id", "gene_name","coding_noncoding_chr")]
genes_candidatesORFs = genes_candidatesORFs %>% subset(coding_noncoding_chr != "other ncORFs") %>% unique()
table(genes_candidatesORFs$coding_noncoding_chr)

temp  = candidatesORFs[,c("chr","coding_noncoding_chr","gene_name")]
temp = temp %>% subset(coding_noncoding_chr != "other ncORFs") %>% subset(grepl("nonX", coding_noncoding_chr)) %>% unique()
table(temp$chr, temp$coding_noncoding_chr)
#############

peaks_annotation_df_candidates <- merge(peaks_annotation_df, genes_candidatesORFs, by.x="feature", by.y="transcript_id", all.x=T, all.y=T)
peaks_annotation_df_candidates <- peaks_annotation_df_candidates %>% unique()          

## only peaks associated to candidates
peaks_annotation_df_candidates <- merge(peaks_annotation_df, genes_candidatesORFs, by.x="feature", by.y="transcript_id")
peaks_annotation_df_candidates <- peaks_annotation_df_candidates %>% unique()     
peaks_annotation_df_candidates$ctype <- gsub("_.*","",peaks_annotation_df_candidates$V4)     
table(peaks_annotation_df_candidates$coding_noncoding_chr)

table(peaks_annotation_df_candidates[,c("gene_name", "coding_noncoding_chr")] %>%  unique() %>% pull("coding_noncoding_chr"))

tp_ctypes <- peaks_annotation_df_candidates[,c("gene_name", "coding_noncoding_chr","ctype")] %>%  unique()
table(tp_ctypes$coding_noncoding_chr,  tp_ctypes$ctype)

## selected cancertypes
selected_tp_ctypes = tp_ctypes %>% subset(ctype %in% c("BRCA","BLCA","LIHC","COAD","KIRC","LUAD","LUSC","PRAD"))
table(selected_tp_ctypes$coding_noncoding_chr,  selected_tp_ctypes$ctype)
