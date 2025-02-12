library(VennDiagram)

pdf(file="/home/marta/Documents/SB_IMMUNOPEP_VEN.pdf")
draw.pairwise.venn(area1 = 10372, ## SB
                   area2 = 254, ## immunopeptidomics
                   cross.area = 114, ## overlap
                   category = c("SB predictions", "immunopeptidomics"),
                   fill = c("#d6512c", "#f9d360"),
                   scaled=T,
                   ext.text=FALSE,
                   print.mode=c("raw"))
dev.off()
