library(Rsamtools)
library(ShortRead)
library(GenomicRanges)


fa = FaFile("/cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fa")

data <- read.csv("/cluster/scratch3/vyp-scratch2/Bochukova_RNASeq/processed/set3_hg38/dexseq/control_PW/Bochukova_set3_control_PW_SignificantExons.csv", stringsAsFactors = FALSE)
data <- data[ !is.na(data$strand), ]


tab.n.exons.5pc <- table(subset(data, FDR < 0.1)$external_gene_id)
tab.n.exons.1pc <- subset(data, FDR < 0.01)$external_gene_id

good.genes <- names(which(tab.n.exons.5pc == 1))  ## good genes have a SINGLE exon with FDR < 5%
v.good.genes <- intersect( good.genes, tab.n.exons.1pc)  ## they also have at least 1 exon with FDR < 1%

## below the table we trust
sig.data <- data[ data$external_gene_id %in% v.good.genes & !is.na(data$external_gene_id) & (data$FDR < 0.01) & !is.na(data$FDR), ]


nseqs <- min(280, nrow(sig.data))
my.width <- 220
my.grange <- GenomicRanges::GRanges(sig.data$chromosome[1:nseqs], IRanges(sig.data$exon.start[1:nseqs]- my.width, sig.data$exon.end[1:nseqs]+my.width))


seq <- getSeq(fa, my.grange)

## now we work with the orientation to make it all consistent
seq.orientated <- ifelse (sig.data$strand[ 1:nseqs ] == 1, as.character(seq), as.character(Biostrings::reverseComplement(seq)))
seq.orientated <- DNAStringSet(seq.orientated)
names(seq.orientated) <- paste0(sig.data$external_gene_id[ 1:nseqs], "_", 1:length(seq.orientated))

seq.orientated.5p <- subseq(seq.orientated, start = 1, end = my.width - 10)
seq.orientated.3p <- subseq(seq.orientated, start = width(seq.orientated) - my.width + 10, end = width(seq.orientated))

writeFasta(seq.orientated.5p, file = "Bochukova_250bp_5p.fa")
writeFasta(seq.orientated.3p, file = "Bochukova_250bp_3p.fa")



print(table(as.character(subseq(seq.orientated, start = my.width-1, end = my.width))))
print(table(as.character(subseq(seq.orientated, start = width(seq.orientated) - my.width + 1, end = width(seq.orientated)- my.width + 2))))
