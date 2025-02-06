# Version D With ChIPpeakAnno
library("ChIPpeakAnno")
library("ChIPseeker")
library("tidyr")
library("dplyr")

gtf_file <- "/data/genomics/marta/genomes/Annot_files_GTF/GENCODE_47_M36/gencode.v47.primary_assembly.annotation.fixed.gtf"
txdb     <- makeTxDbFromGFF(gtf_file)
annotatio_granges <- transcripts(txdb)
names(annotatio_granges) <- mcols(annotatio_granges)$tx_name

atac.peak.files <- c(
  list.files(path="/datasets/marta/scTestis_GSE120508/GSE120507_ATAC", all.files=T, full.names = T,pattern = "narrowPeak.bed")
)


peaks_annotation <- NULL
for (peak_file in atac.peak.files){
  print(peak_file)
  peaks <- readPeakFile(peak_file)
  
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
#############

peaks_annotation_df_candidates <- merge(peaks_annotation_df, genes_candidatesORFs, by.x="feature", by.y="transcript_id", all.x=T, all.y=T)
peaks_annotation_df_candidates <- peaks_annotation_df_candidates %>% unique()          

## only peaks associated to candidates
peaks_annotation_df_candidates <- merge(peaks_annotation_df, genes_candidatesORFs, by.x="feature", by.y="transcript_id")
peaks_annotation_df_candidates <- peaks_annotation_df_candidates %>% unique()          
table(peaks_annotation_df_candidates$coding_noncoding_chr)
