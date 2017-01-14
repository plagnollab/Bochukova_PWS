library(SGSeq) 
library(DEXSeq)

load("data/tx.Rdata") 
load("data/Bochukova_set3_info.RData")
load("data/Bochukova_realigned/sgfc.RData")

support <- read.table("/SAN/biomed/biomed14/vyp-scratch/Bochukova/hg38_kitty/Bochukova_set3.tab", header = T, stringsAsFactor = F) 

load("data/Bochukova_realigned/sgfc_pred.RData") 


## Now take the counts and pretend each event is one gene and each variant is one exon 
## and run in dexseq 

varcounts <- countsVariant(sgvc)
vid <- variantID(sgvc) 
eid <- eventID(sgvc) 



formuladispersion <- ~ sample + (condition + type) * exon
formula0 <-  ~ sample + condition
formula1 <-  ~ sample + condition * exon

DexSeqExons.loc <- DEXSeqDataSet(countData = varcounts,
                                              sampleData = support,
                                              design = formula1, 
					featureID = as.factor(vid),  
					groupID = as.factor(eid)) 

DexSeqExons.loc <- estimateSizeFactors(DexSeqExons.loc)
DexSeqExons.loc <- DEXSeq::estimateDispersions(DexSeqExons.loc)
DexSeqExons.loc <- DEXSeq::testForDEU(DexSeqExons.loc) 
DexSeqExons.loc <- DEXSeq::estimateExonFoldChanges(DexSeqExons.loc) 

save(DexSeqExons.loc, file = "data/Bochukova_realigned/dexseq.RData") 

