library(dplyr)
library(tidyr)

## read table of counts
toc = read.csv("/users/genomics/marta/TestisProject_SaraRazquin/with_TranscriptomeReconstruction/v47/cancers/featureCounts/table_of_counts_TPMs_TCGA-COAD.csv")
toc_long = toc %>% pivot_longer(cols=-c(transcript_id,Length), names_to = "filename", values_to = "TPM")
toc_long$filename = gsub("^X","",toc_long$filename)
toc_long$filename = gsub("\\.","-",toc_long$filename)

## reading patient data
patients = read.csv("/users/genomics/marta/cancers_RNASeq/TCGA_COAD/results/patients.csv")
patients_long = patients %>% pivot_longer(cols=-c(patient), values_to = "filename",names_to = "normal_tumor")
patients_long$patientID = paste0(patients_long$patient,"_",patients_long$normal_tumor)

## substitute filename by patient ID
toc_long_merged = merge(toc_long, patients_long, by="filename")
toc_names = toc_long_merged %>% dplyr::select(transcript_id, Length, TPM, patientID) %>% pivot_wider(names_from = patientID, values_from = "TPM")
