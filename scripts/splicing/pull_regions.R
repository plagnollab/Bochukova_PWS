library(Rsamtools)
library(ShortRead)
library(GenomicRanges)


fa = FaFile("/scratch2/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fa")

data <- read.csv("/scratch2/vyp-scratch2/Bochukova_RNASeq/processed/set3_hg38/dexseq/control_PW/Bochukova_set3_control_PW_SignificantExons.csv")


nseqs <- 300

my.grange <- GenomicRanges::GRanges(data$chromosome[1:nseqs], IRanges(data$exon.start[1:nseqs]-40, data$exon.end[1:nseqs]+40))

my.grange.reduced <- my.grange[  data$exon.end[1:nseqs] - data$exon.start[1:nseqs] < 500 ]

seq = getSeq(fa, my.grange.reduced)
names(seq) <- paste0(names(seq), "_", 1:length(seq))

writeFasta(seq, file = "exons_spliced.fa")
