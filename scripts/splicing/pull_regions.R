library(Rsamtools)
library(ShortRead)
library(GenomicRanges)


fa = FaFile("/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fa")

data <- read.csv("/scratch2/vyp-scratch2/Bochukova_RNASeq/processed/set3_hg38/dexseq/control_PW/Bochukova_set3_control_PW_SignificantExons.csv")
data <- data[ !is.na(data$strand), ]

nseqs <- 300

my.grange <- GenomicRanges::GRanges(data$chromosome[1:nseqs], IRanges(data$exon.start[1:nseqs]-40, data$exon.end[1:nseqs]+40))


seq <- getSeq(fa, my.grange)

seq.orientated <- ifelse (data$strand[ 1:nseqs ] == 1, as.character(seq), as.character(Biostrings::reverseComplement(seq)))
seq.orientated <- DNAStringSet(seq.orientated)

seq.orientated <- seq.orientated[ width(seq.orientated) < 500 ]
names(seq.orientated) <- paste0(names(seq.orientated), "_", 1:length(seq.orientated))

writeFasta(seq.orientated, file = "exons_spliced.fa")


print(table(as.character(subseq(seq, start = 39, end = 40))))
print(table(as.character(subseq(seq.orientated, start = 39, end = 40))))
print(table(as.character(subseq(seq.orientated, start = width(seq.orientated) - 39, end = width(seq.orientated)-38))))
