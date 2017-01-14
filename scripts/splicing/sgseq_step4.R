library(SGSeq) 

## this step is very slow

load("data/tx.Rdata") 
load("data/Bochukova_set3_info.RData")
load("data/Bochukova_realigned/sgfc.RData")

support <- read.table("/SAN/biomed/biomed14/vyp-scratch/Bochukova/hg38_kitty/Bochukova_set3.tab", header = T, stringsAsFactor = F) 

sgvc_pred <- analyzeVariants(sgfc, min_denominator = 10)
sgvc <- getSGVariantCounts(rowRanges(sgvc_pred), sample_info = sample.info, features = sgf ) 
save(sgvc_pred, sgvc, file = "data/Bochukova_realigned/sgfc_pred.RData") 


