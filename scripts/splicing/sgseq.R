library(SGSeq) 

#tx <- importTranscripts("/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/Homo_sapiens.GRCh38.82_fixed.gtf") 
#txf <- convertToTxFeatures(tx) 
#sgf <- convertToSGFeatures(txf) 

#save(tx, txf, sgf, file = "/SAN/biomed/biomed14/vyp-scratch/Bochukova_realigned/tx.Rdata") 
load("/SAN/biomed/biomed14/vyp-scratch/Bochukova_realigned/tx.Rdata") 

#sample.tab <- read.table("/SAN/biomed/biomed14/vyp-scratch/Bochukova_realigned/Bochukova_set3_info.tab", header = T, stringsAsFactor = F) 
support <- read.table("/SAN/biomed/biomed14/vyp-scratch/Bochukova_realigned/Bochukova_set3.tab", header = T, stringsAsFactor = F) 
#sample.info <- getBamInfo(sample.tab) 
#save(sample.info, file = "/SAN/biomed/biomed14/vyp-scratch/Bochukova_realigned/Bochukova_set3_info.RData") 
load("/SAN/biomed/biomed14/vyp-scratch/Bochukova_realigned/Bochukova_set3_info.RData") 

#for testing only 
#chr22 <- GRanges(seqnames = "chr22")
#sgf_chr <- sgf[seqnames(sgf) == "chr22"]
#sgfc <- analyzeFeatures(sample.info, features = sgf)

#save(sgfc, file = "/SAN/biomed/biomed14/vyp-scratch/Bochukova_realigned/sgfc.RData") 
load("/SAN/biomed/biomed14/vyp-scratch/Bochukova_realigned/sgfc.RData") 
#sgvc_pred <- analyzeVariants(sgfc, min_denominator = 10)
#sgvc <- getSGVariantCounts(rowRanges(sgvc_pred), sample_info = sample.info, features = sgf ) 

#save(sgvc_pred, sgvc, file = "/SAN/biomed/biomed14/vyp-scratch/Bochukova_realigned/sgfc_pred.RData") 
load("/SAN/biomed/biomed14/vyp-scratch/Bochukova_realigned/sgfc_pred.RData") 

## Now take the counts and pretend each event is one gene and each variant is one exon 
## and run in dexseq 

varcounts <- countsVariant(sgvc)
vid <- variantID(sgvc) 
eid <- eventID(sgvc) 

library(DEXSeq) 

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

save(DexSeqExons.loc, file = "/SAN/biomed/biomed14/vyp-scratch/Bochukova_realigned/dexseq.RData") 


res <- DEXSeq::DEXSeqResults (DexSeqExons.loc)
res.clean <- as.data.frame(res) 
res.clean$FDR <- p.adjust(res.clean$pvalue, method = 'fdr')
sgvc.df <- as.data.frame(mcols(sgvc) ) 
res.clean$geneName <- sgvc.df$geneName
res.clean$variantType <- sgvc.df$variantType
res.clean$from <- sgvc.df$from
res.clean$to <- sgvc.df$to
res.clean$type <- sgvc.df$type

res.clean <- res.clean[order(res.clean)$pvalue, ] 
write.table(res.clean, file = "sgseq_res.csv", sep = ",") 
