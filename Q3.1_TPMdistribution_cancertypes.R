library(ggplot2)
library(tidyr)
library(dplyr)

cancers = c("BRCA","BLCA","LUAD","KIRC","KIRP","PRAD","LUSC","LIHC")#,"COAD")
tcga_projects=c("TCGA-BRCA","TCGA-LUSC","TCGA-PRAD","TCGA-KIRC","TCGA-KIRP","TCGA-LUAD","TCGA-BLCA")#,"TCGA-LIHC"]
other_projects=c("GSE102101_KIRC","GSE133624_BLCA","GSE22260_PRAD","GSE89223_PRAD","PRJEB2449_PRAD","SRP238334_KIRC","GSE103001_BRCA","GSE214846_LIHC","GSE229705_LUAD","TCGA_COAD")#,"SRP107326_COAD"]
manuscript_projects = c("liver_adjacent_totalRNA_LIHC","hcc_normal_totalRNA_LIHC","GSE193567_LIHC","LIHC_TCGA_LIHC")
all_projects = c(tcga_projects,other_projects,manuscript_projects)

cancer_dir="/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/cancers"
transcript_gene=read.csv("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/human/newReference_Resconstructed/1transcript_1gene.reconstructed.csv")

for (ctype in cancers) {
  print(ctype)
  
  tableofcounts = read.csv(file.path(cancer_dir, paste0("merged_fc_", ctype, ".csv")))
  names(tableofcounts) = gsub("\\.","-", names(tableofcounts))
  bigpatients = read.csv(file.path(cancer_dir, paste0("merged_patients_", ctype, ".csv")))
  bigpatients = bigpatients %>% select(tumor, project)
  names(bigpatients) = c("sample","project")
  
  tableofcounts_tumor = tableofcounts %>% select(transcript_id, Length, all_of(bigpatients$sample))
  ## add gene type data
  tableofcounts_tumor = merge(tableofcounts_tumor, transcript_gene, by="transcript_id")
  ## subset interesting gene types
  tableofcounts_tumor = tableofcounts_tumor %>% subset(gene_type %in% c("lncRNA","novel","processed_pseudogene","protein_coding"))
  tableofcounts_tumor_long = tableofcounts_tumor %>% pivot_longer(cols=-c(transcript_id, gene_id, gene_name, gene_type, Length), names_to = "sample", values_to = "TPM") 
  ## add project data
  combined = merge(tableofcounts_tumor_long, bigpatients, by="sample")
  combined$logTPM = log10(combined$TPM)
  
  ggplot(combined, aes(x=logTPM, color=project)) +
    geom_density() +
    geom_vline(xintercept = 1) +
    labs(title=ctype,
         x="log10TPM") +
    theme_classic() +
    facet_wrap(~ gene_type)
  ggsave(paste0("/projects_eg/projects/marta/TestisRestricted_Microproteins_TSA/with_TranscriptomeReconstruction/plots/PNG/TPMdistribution_",ctype,".png"), width=5.41, height=4.66)
  ggsave(paste0("/projects_eg/projects/marta/TestisRestricted_Microproteins_TSA/with_TranscriptomeReconstruction/plots/PDF/TPMdistribution_",ctype,".pdf"), width=5.41, height=4.66)
  
}
