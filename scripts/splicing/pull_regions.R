library(Rsamtools)
library(ShortRead)
library(GenomicRanges)


fa = FaFile("/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fa")

#data <- read.csv("/scratch2/vyp-scratch2/Bochukova_RNASeq/processed/set3_hg38/dexseq/control_PW/Bochukova_set3_control_PW_SignificantExons.csv")
#data <- data[ !is.na(data$strand), ]

nseqs <- 280
my.width <- 220
my.grange <- GenomicRanges::GRanges(data$chromosome[1:nseqs], IRanges(data$exon.start[1:nseqs]- my.width, data$exon.end[1:nseqs]+my.width))


seq <- getSeq(fa, my.grange)

seq.orientated <- ifelse (data$strand[ 1:nseqs ] == 1, as.character(seq), as.character(Biostrings::reverseComplement(seq)))
seq.orientated <- DNAStringSet(seq.orientated)
names(seq.orientated) <- paste0(data$external_gene_id[ 1:nseqs], "_", 1:length(seq.orientated))

seq.orientated.5p <- subseq(seq.orientated, start = 1, end = my.width - 10)
seq.orientated.3p <- subseq(seq.orientated, start = width(seq.orientated) - my.width + 10, end = width(seq.orientated))

writeFasta(seq.orientated.5p, file = "Bochukova_250bp_5p.fa")
writeFasta(seq.orientated.3p, file = "Bochukova_250bp_3p.fa")



print(table(as.character(subseq(seq.orientated, start = my.width-1, end = my.width))))
print(table(as.character(subseq(seq.orientated, start = width(seq.orientated) - my.width + 1, end = width(seq.orientated)- my.width + 2))))
