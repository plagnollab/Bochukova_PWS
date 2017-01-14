library(DESeq2)

targets <- read.csv("data/SNORD116_targets.csv", stringsAsFactors = FALSE)
exons.sig.v2 <- read.csv("sgseq_res.csv", stringsAsFactors = FALSE)


all.genes <- data.frame(gene = unique(exons.sig.v2$external_geneName))
all.genes$signif <- all.genes$gene %in% exons.sig.v2$external_geneName[ exons.sig.v2$FDR < 0.1 ]
all.genes$SNO.target <- all.genes$gene %in% targets$Gene.name
fisher.test(table(all.genes$signif, all.genes$SNO.target))



stop()
load("set3_hg38/dexseq/control_PW/dexseq_Bochukova_set3_control_PW.RData")
exons.sig <- read.csv("set3_hg38/dexseq/control_PW/Bochukova_set3_control_PW_SignificantExons.csv")
exons.sig <- exons.sig[ ! grepl(pattern = "^RP|^AC|^AL", exons.sig$external_gene_id), ]
exons.sig <- exons.sig[ !is.na(exons.sig$external_gene_id), ]
exons.sig <- exons.sig[ exons.sig$meanBase > 1, ]



all.genes <- data.frame(gene = unique(exons.sig$external_gene_id))
all.genes$signif <- all.genes$gene %in% exons.sig$external_gene_id[ exons.sig$FDR < 0.01 ]
all.genes$SNO.target <- all.genes$gene %in% targets$Gene.name

write.csv(x = all.genes, row.names = FALSE, file = "splicing_SNORD116.csv")

print(fisher.test(table(all.genes$signif, all.genes$SNO.target)))
