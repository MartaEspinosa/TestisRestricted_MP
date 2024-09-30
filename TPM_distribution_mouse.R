library(dplyr)
library(tidyr)
library(ggplot2)

# wd = "/users/genomics/marta/TestisProject_SaraRazquin"
wd = "~/Downloads/IMIM"

biomart = read.csv(file.path(wd,"with_TranscriptomeReconstruction/mouse/newReference_Resconstructed/transcript_gene.csv"))
biomart$gene_id = gsub("\\..*" biomart$gene_id)
biomart$transcript_id = gsub("\\..*" biomart$transcript_id)

TPMs = read.csv(file.path(wd,"with_TranscriptomeReconstruction/mouse/featureCounts_gffcompare/table_of_counts_TPMs.csv"))
names(TPMs)[1] = "transcript_id"

TPMs_complete = merge(TPMs, biomart, by="transcript_id")
TPMs_complete = TPMs_complete %>% pivot_longer(cols=-c(gene_id, transcript_id, gene_name, gene_type), names_to = "sample values_to = "TPM")
TPMs_complete$logTPM = log(TPMs_complete$TPM)

ggplot(TPMs_complete %>% subset(gene_type == "lncRNA" | gene_type == "processed_pseudogene" | gene_type == "novel" | gene_type == "protein_coding"), aes(x=logTPM, fill=gene_type)) +
  geom_density(alpha=.5) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  theme(legend.position = "top") +
  facet_wrap(~ sample, ncol=3)
ggsave(file.path(wd,"with_TranscriptomeReconstruction/mouse/plots/PNG/TPM_density.png"), width=6.64, height=8.78)
ggsave(file.path(wd,"with_TranscriptomeReconstruction/mouse/plots/PDF/TPM_density.pdf"), width=6.64, height=8.78)

ggplot(TPMs_complete %>% subset(gene_type == "lncRNA" | gene_type == "processed_pseudogene" | gene_type == "novel" | gene_type == "protein_coding"), aes(x=gene_type, y=logTPM, fill=gene_type)) +
  geom_boxplot() +
  geom_vline(xintercept = 0) +
  theme_classic() +
  theme(legend.position = "top
        axis.text.x = element_text(angle=45, vjust=.75)) +
  facet_wrap(~ sample, ncol=3)
ggsave(file.path(wd,"with_TranscriptomeReconstruction/mouse/plots/PNG/TPM_boxplot.png"), width=6.64, height=8.78)
ggsave(file.path(wd,"with_TranscriptomeReconstruction/mouse/plots/PDF/TPM_boxplot.pdf"), width=6.64, height=8.78)

## length distribution
toc = read.csv(file.path(wd,"with_TranscriptomeReconstruction/mouse/featureCounts_gffcompare/gffcompare_stranded_featureCounts.txt"), sep="\t skip = 1)
names(toc)[1] = "transcript_id"
toc$transcript_id = gsub("\\..*"toc$transcript_id)
toc = toc %>% select(transcript_id, Chr, Length)
TPMs_length = merge(TPMs_complete, toc, by="transcript_id")

ggplot(TPMs_length %>% subset(gene_type == "lncRNA" | gene_type == "processed_pseudogene" | gene_type == "novel" | gene_type == "protein_coding"), aes(x=Length, fill=gene_type)) +
  geom_density(alpha=.5) +
  scale_x_continuous(trans="log10") +
  geom_vline(xintercept=300) +
  theme_classic() +
  theme(legend.position = "top") +
  facet_wrap(~ sample, ncol=3)
ggsave(file.path(wd,"with_TranscriptomeReconstruction/mouse/plots/PNG/length_tableofcounts_density.png"), height=6.64, width=8.78)
ggsave(file.path(wd,"with_TranscriptomeReconstruction/mouse/plots/PDF/length_tableofcounts_density.pDf"), height=6.64, width=8.78)

TPMs_length_300 = TPMs_length %>%
  filter(!(gene_type == "novel" & Length < 300))

ggplot(TPMs_length_300 %>% subset(gene_type == "lncRNA" | gene_type == "processed_pseudogene" | gene_type == "novel" | gene_type == "protein_coding"), aes(x=Length, fill=gene_type)) +
  geom_density(alpha=.5) +
  scale_x_continuous(trans="log10") +
  geom_vline(xintercept=300) +
  theme_classic() +
  theme(legend.position = "top") +
  facet_wrap(~ sample, ncol=3)
ggsave(file.path(wd,"with_TranscriptomeReconstruction/mouse/plots/PNG/length_tableofcounts_density.300.png"), height=6.64, width=8.78)
ggsave(file.path(wd,"with_TranscriptomeReconstruction/mouse/plots/PDF/length_tableofcounts_density.300.pdf"), height=6.64, width=8.78)

ggplot(TPMs_length_300 %>% subset(gene_type == "lncRNA" | gene_type == "processed_pseudogene" | gene_type == "novel" | gene_type == "protein_coding"), aes(x=logTPM, fill=gene_type)) +
  geom_density(alpha=.5) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  theme(legend.position = "top") +
  facet_wrap(~ sample, ncol=3)
ggsave(file.path(wd,"with_TranscriptomeReconstruction/mouse/plots/PNG/TPM_density.300.png"), width=6.64, height=8.78)
ggsave(file.path(wd,"with_TranscriptomeReconstruction/mouse/plots/PDF/TPM_density.300.pdf"), width=6.64, height=8.78)

ggplot(TPMs_length_300 %>% subset(gene_type == "lncRNA" | gene_type == "processed_pseudogene" | gene_type == "novel" | gene_type == "protein_coding"), aes(x=gene_type, y=logTPM, fill=gene_type)) +
  geom_boxplot() +
  geom_vline(xintercept = 0) +
  theme_classic() +
  theme(legend.position = "top
        axis.text.x = element_text(angle=45, vjust=.75)) +
  facet_wrap(~ sample, ncol=3)
ggsave(file.path(wd,"with_TranscriptomeReconstruction/mouse/plots/PNG/TPM_boxplot.300.png"), width=6.64, height=8.78)
ggsave(file.path(wd,"with_TranscriptomeReconstruction/mouse/plots/PDF/TPM_boxplot.300.pdf"), width=6.64, height=8.78)

### shared?
shared = TPMs_length_300 %>% 
  subset(gene_type == "lncRNA" | gene_type == "processed_pseudogene" | gene_type == "novel" | gene_type == "protein_coding") %>%
  subset(TPM > 1) %>%
  group_by(gene_id, gene_type) %>%
  summarise(num_patients=n())

shared_n = shared %>% group_by(gene_type, num_patients) %>% summarise(num_genes = n())

ggplot(shared_n, aes(x=num_patients, y=num_genes, color=gene_type)) +
  geom_point(stat="identity") +
  geom_line()+
  scale_y_continuous(trans="log10") +
  ggtitle("Shared genes, all samples | > 300 pb & > 1 TPM") +
  theme_classic()
ggsave(file.path(wd,"with_TranscriptomeReconstruction/mouse/plots/PNG/shared_genes.300.png"), width=6.64, height=8.78)
ggsave(file.path(wd,"with_TranscriptomeReconstruction/mouse/plots/PDF/shared_genes.300.pdf"), width=6.64, height=8.78)

