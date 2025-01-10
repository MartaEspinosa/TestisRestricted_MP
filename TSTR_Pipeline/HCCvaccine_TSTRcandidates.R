## Candidates from testis for HCC-vaccine

library(dplyr)
library(tidyr)

## setup
tumorReactDIR = "/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/cancers/log2ratio3x/cancertypes"
annot = read.csv("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/human/newReference_Resconstructed/transID_geneID_isoforms_selected.1to1.csv")


## SciAdv
SciAdv = read.csv("/projects_eg/projects/marta/table_to_heatmap_noabundantcase3_IEAtlas.csv")
SciAdv = SciAdv %>% select(gene_id, n, RiboSeq) %>% unique()
names(SciAdv)[2] = "n_SciAdv"

## log2ratio 
HCC_log2ratio = read.csv("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/cancers/log2ratio3x/cancertypes/LIHC/log2ratio3x.csv")
HCC_log2ratio = HCC_log2ratio %>% select(gene_id, log2ratio)

## other immuno
# by peptide
immuno = read.csv("/projects_eg/projects/marta/TestisRestricted_Microproteins_TSA/with_TranscriptomeReconstruction/Multimap_altORFs/Q5_immunopeptidomics/human/TOv3x_5percent_TestisRestrictedGTEx_Translated_Ctypes_log2ratio3x_immuno.csv")
immuno = immuno %>% select(Peptide, gene_id_x, geneORFtype, source, matching_ORFpep)
names(immuno) = c("Peptide_broadSources","gene_id","geneORFtype", "BroadSource","immuno_ORFpep")

# by gene
# immuno = read.csv("/projects_eg/projects/marta/TestisRestricted_Microproteins_TSA/with_TranscriptomeReconstruction/Multimap_altORFs/Q5_immunopeptidomics/human/TOv3x_5percent_TestisRestrictedGTEx_Translated_Ctypes_log2ratio3x_immunobyGene.csv")
# immuno = immuno %>% select(Peptide, gene_id, geneORFtype, source)
# names(immuno) = c("Peptide_broadSources","gene_id","geneORFtype", "BroadSource")

# ## genes
# candidates = read.csv(file.path(tumorReactDIR,"TOv3x_5percent_TestisRestrictedGTEx_Translated_Ctypes_log2ratio3xMEANgenes.csv"))
# candidates = merge(candidates, annot %>% select(gene_id, chr), by="gene_id")
# candidates = candidates %>% mutate(coding_noncoding_chr = case_when(gene_type == "protein_coding" & chr == "X" ~ "CT-X",
#                                                                     gene_type == "protein_coding" & chr != "X" ~ "CT-nonX",
#                                                                     gene_type != "protein_coding" & chr == "X" ~ "Noncoding-X",
#                                                                     gene_type != "protein_coding" & chr != "X" ~ "Noncoding-nonX"),
#                                    type = case_when(gene_type == "protein_coding" ~ "CT",
#                                                     TRUE ~ "Noncoding"))

## ORFs
candidatesORFs = read.csv(file.path(tumorReactDIR,"TOv3x_5percent_TestisRestrictedGTEx_Translated_Ctypes_log2ratio3xMEAN.csv"))
candidatesORFs = merge(candidatesORFs, annot %>% select(gene_id, chr), by="gene_id")
candidatesORFs = candidatesORFs %>% mutate(coding_noncoding_chr = case_when(gene_type == "protein_coding" & chr == "X" ~ "CT-X",
                                                                        gene_type == "protein_coding" & chr != "X" ~ "CT-nonX",
                                                                        gene_type != "protein_coding" & chr == "X" ~ "Noncoding-X",
                                                                        gene_type != "protein_coding" & chr != "X" ~ "Noncoding-nonX"),
                                           type = case_when(gene_type == "protein_coding" ~ "CT",
                                                            TRUE ~ "Noncoding"))

## Expression of gene in tumor
TSTR_tumor = read.csv(file.path(tumorReactDIR, "TSTR_Expression_tumor.csv"))
stats_TSTR_tumor = TSTR_tumor %>% group_by(gene_id) %>% summarise(medianTPM_tumor = median(TPM), maxTPM_tumor = max(TPM))
TSTR_normal = read.csv(file.path(tumorReactDIR, "TSTR_Expression_normal.csv"))
stats_TSTR_normal = TSTR_normal %>% group_by(gene_id) %>% summarise(medianTPM_normal = median(TPM), maxTPM_normal = max(TPM))

stats_TSTR = merge(stats_TSTR_normal, stats_TSTR_tumor, by="gene_id")

######### HCC - selection ######### 
candidatesORFs_HCC = candidatesORFs %>% subset(ctype == "LIHC")
candidatesORFs_HCC = candidatesORFs_HCC %>% subset(geneORFtype != "protein_coding_canonical")
table(candidatesORFs_HCC$geneORFtype, candidatesORFs_HCC$coding_noncoding_chr)

candidatesORFs_pancancer = candidatesORFs %>% select(gene_id, ctype) %>% unique() %>% subset(gene_id %in% candidatesORFs_HCC$gene_id)
candidatesORFs_pancancer = candidatesORFs_pancancer %>% group_by(gene_id) %>% count()
names(candidatesORFs_pancancer)[2] = "num_ctypes"

candORFs_HCC = merge(candidatesORFs_HCC, candidatesORFs_pancancer, on="gene_id")
candORFs_HCC_SciAdv = merge(candORFs_HCC, SciAdv, on="gene_id", all.x=T)
candORFs_HCC_SciAdv = merge(candORFs_HCC_SciAdv, HCC_log2ratio, on="gene_id", all.x=T)
# candORFs_HCC_SciAdv = merge(candORFs_HCC_SciAdv, immuno, on=c("gene_id","ORFpep"), all.x=T)
candORFs_HCC_SciAdv = merge(candORFs_HCC_SciAdv, immuno, on=c("gene_id","geneORFtype"), all.x=T)
candORFs_HCC_SciAdv = merge(candORFs_HCC_SciAdv, stats_TSTR, on="gene_id", all.x=T)

## add immunopeptidomics
########### BEIJER ###########
Beijer_5percent = read.csv("/users/genomics/marta/241205_Immunopeptidomics_HCC/HCC/FDR5percent/FDR5percent_TSA_TSTR_analysis.csv")
Beijer_5percent$source = "Beijer_5percent"
Beijer_5percent = Beijer_5percent %>% select(sequence, sample, accessions, gene_id, source)

########### BEIJER - PT ###########
BeijerPT_5percent = read.csv("/users/genomics/marta/241205_Immunopeptidomics_HCC/HCC_PT/FDR5percent/FDR5percent_TSA_TSTR_analysis.csv")
BeijerPT_5percent$source = "BeijerPT_5percent"
BeijerPT_5percent = BeijerPT_5percent %>% select(sequence, sample, accessions, gene_id, source)

########### cell lines ###########
celllines_5percent = read.csv("/users/genomics/marta/241205_Immunopeptidomics_HCC/celllines/FDR5percent/FDR5percent_TSA_TSTR_analysis.csv")
celllines_5percent$source = "celllines_5percent"
celllines_5percent = celllines_5percent %>% select(sequence, sample, accessions, gene_id, source)

########### Löffer ###########
Lof_5percent = read.csv("/users/genomics/marta/241205_Immunopeptidomics_HCC/Löffer/FDR5percent/FDR5percent_TSA_TSTR_analysis.csv")
Lof_5percent$source = "Löffer_5percent"
Lof_5percent = Lof_5percent %>% select(sequence, sample, accessions, gene_id, source)

MHCquant = rbind(celllines_5percent, BeijerPT_5percent, Beijer_5percent, Lof_5percent)

candORFs_HCC_SciAdv_MHCq = merge(candORFs_HCC_SciAdv, MHCquant, all.x=T, on="gene_id")
candORFs_HCC_SciAdv_MHCq = candORFs_HCC_SciAdv_MHCq %>% select(gene_id, gene_name, geneORFtype, percentage_num_patients_overexpr, n_SciAdv, log2ratio, coding_noncoding_chr, source, sequence, BroadSource, Peptide_broadSources, orfID, accessions, ORFpep, everything()) %>% unique()
write.csv(candORFs_HCC_SciAdv_MHCq, "/home/marta/candORFs_HCC_SciAdv_MHCq.csv", row.names = F)

### lncRNAs
candidatesORFs_HCC_lncRNA = candORFs_HCC_SciAdv_MHCq %>% subset(grepl("lncRNA", geneORFtype))

### novel
candidatesORFs_HCC_novel = candORFs_HCC_SciAdv_MHCq %>% subset(grepl("novel", geneORFtype))

### altORFs
candidatesORFs_HCC_altORFs = candORFs_HCC_SciAdv_MHCq %>% subset(grepl("ORF", geneORFtype))

