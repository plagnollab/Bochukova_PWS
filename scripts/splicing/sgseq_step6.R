library(SGSeq) 


support <- read.table("/SAN/biomed/biomed14/vyp-scratch/Bochukova/hg38_kitty/Bochukova_set3.tab", header = T, stringsAsFactor = F) 
load("data/Bochukova_realigned/sgfc_pred.RData", verbose = TRUE)
load("data/Bochukova_realigned/dexseq.RData", verbose = TRUE)
load("data/Bochukova_realigned/dexseq.RData")


res <- DEXSeq::DEXSeqResults (DexSeqExons.loc)
res.clean <- as.data.frame(res)

sgvc.df <- as.data.frame(mcols(sgvc) ) 


res.clean$geneName <- sgvc.df$geneName
res.clean$variantType <- sgvc.df$variantType
res.clean$from <- sgvc.df$from
res.clean$to <- sgvc.df$to
res.clean$type <- sgvc.df$type

res.clean <- res.clean[ !is.na(res.clean$log2fold_PW_control), ]
res.clean$FDR <- p.adjust(res.clean$pvalue, method = 'fdr')


res.clean <- res.clean[order(res.clean$pvalue), ]
res.clean$geneName <- sapply(res.clean$geneName, paste, collapse = "_")

res.clean$external_geneName <- annotation$external_gene_name[ match(res.clean$geneName, table = annotation$EnsemblID) ]

res.clean$variantType.clean <- sapply(res.clean$variantType, paste, collapse = "_")

annotation <- read.table("/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/human_hg38/biomart/biomart_annotations_human.tab", stringsAsFactors = FALSE, header = TRUE)

write.csv(res.clean, file = "sgseq_res.csv",  row.names = FALSE) 


countEvents <- sort(table( dplyr::filter(res.clean, FDR < 0.1)$variantType.clean))

