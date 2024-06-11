library(dplyr)
library(tidyr)
library(ggplot2)


biomart = read.csv("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/macaca/newReference_Resconstructed/transcript_gene.csv")
TPMs = read.csv("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/macaca/featureCounts_gffcompare/table_of_counts_TPMs.csv")
names(TPMs)[1] = "transcript_id"

TPMs_complete = merge(TPMs, biomart, by="transcript_id")
TPMs_complete = TPMs_complete %>% pivot_longer(cols=-c(gene_id, transcript_id, gene_name, gene_type), names_to = "sample", values_to = "TPM")
TPMs_complete$logTPM = log(TPMs_complete$TPM)

ggplot(TPMs_complete %>% subset(gene_type == "lncRNA" | gene_type == "processed_pseudogene" | gene_type == "novel" | gene_type == "protein_coding"), aes(x=logTPM, fill=gene_type)) +
  geom_density(alpha=.5) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  theme(legend.position = "top") +
  facet_wrap(~ sample, ncol=3)
ggsave("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/macaca/plots/PNG/TPM_density.png", height=6.64, width=8.78)
ggsave("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/macaca/plots/PDF/TPM_density.pdf", height=6.64, width=8.78)

ggplot(TPMs_complete %>% subset(gene_type == "lncRNA" | gene_type == "processed_pseudogene" | gene_type == "novel" | gene_type == "protein_coding"), aes(x=gene_type, y=logTPM, fill=gene_type)) +
  geom_boxplot() +
  geom_vline(xintercept = 0) +
  theme_classic() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle=45, vjust=.75)) +
  facet_wrap(~ sample, ncol=3)
ggsave("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/macaca/plots/PNG/TPM_boxplot.png", height=6.64, width=8.78)
ggsave("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/macaca/plots/PDF/TPM_boxplot.pdf", height=6.64, width=8.78)

## length distribution
toc = read.csv("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/macaca/featureCounts_gffcompare/gffcompare_stranded_featureCounts.txt", sep="\t", skip = 1)
names(toc)[1] = "transcript_id"
toc = toc %>% select(transcript_id, Length)
TPMs_length = merge(TPMs_complete, toc, by="transcript_id")

ggplot(TPMs_length %>% subset(gene_type == "lncRNA" | gene_type == "processed_pseudogene" | gene_type == "novel" | gene_type == "protein_coding"), aes(x=Length, fill=gene_type)) +
  geom_density(alpha=.5) +
  scale_x_continuous(trans="log10") +
  theme_classic() +
  theme(legend.position = "top") +
  facet_wrap(~ sample, ncol=3)
ggsave("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/macaca/plots/PNG/length_tableofcounts_density.png", height=6.64, width=8.78)
ggsave("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/macaca/plots/PDF/length_tableofcounts_density.pdf", height=6.64, width=8.78)
