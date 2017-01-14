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
write.table(res.clean, file = "sgseq_res.csv", sep = ",") 
