






library(SGSeq) 


load("data/tx.Rdata", verbose = TRUE)
load("data/Bochukova_set3_info.RData", verbose = TRUE)

sgfc <- analyzeFeatures(sample.info, features = sgf)
save(sgfc, file = "data/Bochukova_realigned/sgfc.RData") 
