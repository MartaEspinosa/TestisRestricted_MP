## Candidates from testis for HCC-vaccine

library(dplyr)
library(tidyr)

## setup
tumorReactDIR = "/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/cancers/log2ratio3x/cancertypes"
annot = read.csv("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/human/newReference_Resconstructed/transID_geneID_isoforms_selected.1to1.csv")


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
TSTR_normal = read.csv(file.path(tumorReactDIR, "TSTR_Expression_normal.csv"))

######### HCC - selection ######### 
candidatesORFs_HCC = candidatesORFs %>% subset(ctype == "LIHC")
candidatesORFs_HCC = candidatesORFs_HCC %>% subset(geneORFtype != "protein_coding_canonical")
table(candidatesORFs_HCC$geneORFtype, candidatesORFs_HCC$coding_noncoding_chr)

candidatesORFs_pancancer = candidatesORFs %>% select(gene_id, ctype) %>% subset(gene_id %in% candidatesORFs_HCC$gene_id)
candidatesORFs_pancancer = candidatesORFs_pancancer %>% group_by(gene_id) %>% count()
names(candidatesORFs_pancancer)[2] = "num_ctypes"

candidatesORFs_HCC = merge(candidatesORFs_HCC, candidatesORFs_pancancer, on="gene_id")
