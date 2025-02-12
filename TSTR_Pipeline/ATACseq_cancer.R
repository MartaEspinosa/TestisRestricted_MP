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
  
  # ann.tmp <- annotatePeakInBatch(myPeakList = peaks,AnnotationData = annotatio_granges ,output = "nearestBiDirectionalPromoters",multiple = T, bindingRegion = c(-1000, 100))
  ann.tmp <- annotatePeakInBatch(myPeakList = peaks,AnnotationData = annotatio_granges ,output = "upstream&inside",multiple = T)
  
  peaks_annotation  <- c(peaks_annotation , ann.tmp)
}

peaks_annotation <- do.call(c, peaks_annotation)
peaks_annotation_df <- as.data.frame(peaks_annotation, row.names = NULL)
peaks_annotation_df$feature <- gsub("\\..","",peaks_annotation_df$feature)
write.csv(peaks_annotation_df, "/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/ATACseq_cancer/AnnotatedPeaks_Upstream.csv")


### Candidates 
tumorReactDIR = "/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/cancers/log2ratio3x/cancertypes/onlyStep1"
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
peaks_annotation_df_candidates$ctype_ATAC <- gsub("_.*","",peaks_annotation_df_candidates$V4)     
table(peaks_annotation_df_candidates$coding_noncoding_chr)
write.csv(peaks_annotation_df_candidates, "/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/ATACseq_cancer/CandidatesORFs_AnnotatedPeaks_Upstream.csv")

table(peaks_annotation_df_candidates[,c("gene_name", "coding_noncoding_chr")] %>%  unique() %>% pull("coding_noncoding_chr"))

tp_ctypes <- peaks_annotation_df_candidates[,c("gene_name", "coding_noncoding_chr","ctype_ATAC")] %>%  unique()
table(tp_ctypes$coding_noncoding_chr,  tp_ctypes$ctype_ATAC)

## same ATAC and expression ctype
selected_tp_ctypes = tp_ctypes %>% subset(ctype_ATAC %in% c("BRCA","BLCA","LIHC","COAD","KIRC","LUAD","LUSC","PRAD"))
table(selected_tp_ctypes$coding_noncoding_chr,  selected_tp_ctypes$ctype_ATAC)

selected_tp_ctypes_same = merge(selected_tp_ctypes, candidatesORFs[,c("gene_name", "ctype")], by.x=c("gene_name","ctype_ATAC"), by.y=c("gene_name","ctype"))
selected_tp_ctypes_same = selected_tp_ctypes_same %>% unique()
table(selected_tp_ctypes_same$coding_noncoding_chr,  selected_tp_ctypes_same$ctype_ATAC)


peaks_annotation_df_candidates_cancertypes_matched = merge(peaks_annotation_df_candidates, selected_tp_ctypes_same, by=c("ctype_ATAC","gene_name", "coding_noncoding_chr"))
write.csv(peaks_annotation_df_candidates_cancertypes_matched, "/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/ATACseq_cancer/CandidatesORFsSameCancerType_AnnotatedPeaks_Upstream.csv")

### GWAS HCC data
# gwas_hcc = read.csv("/datasets/marta/GWAS/j_liver_xintra_cancer.tsv", sep="\t")
# hcc_candidates = candidatesORFs %>% subset(ctype == "LIHC")
# peaks_candidates_hcc = peaks_annotation_df_candidates %>% subset(gene_name %in% hcc_candidates$gene_name) %>% subset(ctype_ATAC == "LIHC")
# 
# snps_in_peaks <- gwas_hcc[,c("chrCHR","POS","SNP","FREQ_European")] %>%
#   inner_join(peaks_candidates_hcc, by = c("chrCHR" = "seqnames")) %>%
#   filter(POS >= start & POS <= end)
