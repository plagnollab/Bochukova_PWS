library(GenomicRanges)
library(Rsamtools)
library(Biostrings)

editing.sites <- read.table("data/Human_AG_NonRepetitive_hg19_v2.txt", header = TRUE, stringsAsFactors = FALSE)
editing.sites <- subset(editing.sites, !annot1 %in% c("intergenic", "intronic") & chromosome %in% paste0("chr", as.character(1:22)))
editing.sites$chromosome <- gsub(pattern = "chr", editing.sites$chromosome, replacement = "")
editing.sites <- editing.sites[ order(as.numeric(editing.sites$chromosome), editing.sites$position), ]


my.ranges <- GRanges(gsub(pattern = "^chr", replacement = "", editing.sites$chromosome), IRanges(editing.sites$position, editing.sites$position ))


fasta.file <- "~/vyp-scratch2/reference_datasets/human_reference_sequence/human_g1k_v37.fasta"
fa = Rsamtools::FaFile(fasta.file)
refSeq <- Rsamtools::getSeq(fa, my.ranges)


print(table(as.character(refSeq), editing.sites$strand))

my.vcf <- data.frame(chromosome = editing.sites$chromosome,
                     position = editing.sites$position,
                     ID = ".",
                     REF = ifelse (editing.sites$strand == "+", "A", "T"),
                     ALT = ifelse (editing.sites$strand == "+", "G", "C"),
                     QUAL = 999,
                     FILTER = "PASS",
                     INFO = ".",
                     stringsAsFactors = FALSE)

cat("##fileformat=VCFv4.1
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=2,length=243199373,assembly=b37>
##contig=<ID=3,length=198022430,assembly=b37>
##contig=<ID=4,length=191154276,assembly=b37>
##contig=<ID=5,length=180915260,assembly=b37>
##contig=<ID=6,length=171115067,assembly=b37>
##contig=<ID=7,length=159138663,assembly=b37>
##contig=<ID=8,length=146364022,assembly=b37>
##contig=<ID=9,length=141213431,assembly=b37>
##contig=<ID=10,length=135534747,assembly=b37>
##contig=<ID=11,length=135006516,assembly=b37>
##contig=<ID=12,length=133851895,assembly=b37>
##contig=<ID=13,length=115169878,assembly=b37>
##contig=<ID=14,length=107349540,assembly=b37>
##contig=<ID=15,length=102531392,assembly=b37>
##contig=<ID=16,length=90354753,assembly=b37>
##contig=<ID=17,length=81195210,assembly=b37>
##contig=<ID=18,length=78077248,assembly=b37>
##contig=<ID=19,length=59128983,assembly=b37>
##contig=<ID=20,length=63025520,assembly=b37>
##contig=<ID=21,length=48129895,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
##contig=<ID=X,length=155270560,assembly=b37>
##contig=<ID=Y,length=59373566,assembly=b37>
##contig=<ID=MT,length=16569,assembly=b37>
##contig=<ID=GL000207.1,length=4262,assembly=b37>
##contig=<ID=GL000226.1,length=15008,assembly=b37>
##contig=<ID=GL000229.1,length=19913,assembly=b37>
##contig=<ID=GL000231.1,length=27386,assembly=b37>
##contig=<ID=GL000210.1,length=27682,assembly=b37>
##contig=<ID=GL000239.1,length=33824,assembly=b37>
##contig=<ID=GL000235.1,length=34474,assembly=b37>
##contig=<ID=GL000201.1,length=36148,assembly=b37>
##contig=<ID=GL000247.1,length=36422,assembly=b37>
##contig=<ID=GL000245.1,length=36651,assembly=b37>
##contig=<ID=GL000197.1,length=37175,assembly=b37>
##contig=<ID=GL000203.1,length=37498,assembly=b37>
##contig=<ID=GL000246.1,length=38154,assembly=b37>
##contig=<ID=GL000249.1,length=38502,assembly=b37>
##contig=<ID=GL000196.1,length=38914,assembly=b37>
##contig=<ID=GL000248.1,length=39786,assembly=b37>
##contig=<ID=GL000244.1,length=39929,assembly=b37>
##contig=<ID=GL000238.1,length=39939,assembly=b37>
##contig=<ID=GL000202.1,length=40103,assembly=b37>
##contig=<ID=GL000234.1,length=40531,assembly=b37>
##contig=<ID=GL000232.1,length=40652,assembly=b37>
##contig=<ID=GL000206.1,length=41001,assembly=b37>
##contig=<ID=GL000240.1,length=41933,assembly=b37>
##contig=<ID=GL000236.1,length=41934,assembly=b37>
##contig=<ID=GL000241.1,length=42152,assembly=b37>
##contig=<ID=GL000243.1,length=43341,assembly=b37>
##contig=<ID=GL000242.1,length=43523,assembly=b37>
##contig=<ID=GL000230.1,length=43691,assembly=b37>
##contig=<ID=GL000237.1,length=45867,assembly=b37>
##contig=<ID=GL000233.1,length=45941,assembly=b37>
##contig=<ID=GL000204.1,length=81310,assembly=b37>
##contig=<ID=GL000198.1,length=90085,assembly=b37>
##contig=<ID=GL000208.1,length=92689,assembly=b37>
##contig=<ID=GL000191.1,length=106433,assembly=b37>
##contig=<ID=GL000227.1,length=128374,assembly=b37>
##contig=<ID=GL000228.1,length=129120,assembly=b37>
##contig=<ID=GL000214.1,length=137718,assembly=b37>
##contig=<ID=GL000221.1,length=155397,assembly=b37>
##contig=<ID=GL000209.1,length=159169,assembly=b37>
##contig=<ID=GL000218.1,length=161147,assembly=b37>
##contig=<ID=GL000220.1,length=161802,assembly=b37>
##contig=<ID=GL000213.1,length=164239,assembly=b37>
##contig=<ID=GL000211.1,length=166566,assembly=b37>
##contig=<ID=GL000199.1,length=169874,assembly=b37>
##contig=<ID=GL000217.1,length=172149,assembly=b37>
##contig=<ID=GL000216.1,length=172294,assembly=b37>
##contig=<ID=GL000215.1,length=172545,assembly=b37>
##contig=<ID=GL000205.1,length=174588,assembly=b37>
##contig=<ID=GL000219.1,length=179198,assembly=b37>
##contig=<ID=GL000224.1,length=179693,assembly=b37>
##contig=<ID=GL000223.1,length=180455,assembly=b37>
##contig=<ID=GL000195.1,length=182896,assembly=b37>
##contig=<ID=GL000212.1,length=186858,assembly=b37>
##contig=<ID=GL000222.1,length=186861,assembly=b37>
##contig=<ID=GL000200.1,length=187035,assembly=b37>
##contig=<ID=GL000193.1,length=189789,assembly=b37>
##contig=<ID=GL000194.1,length=191469,assembly=b37>
##contig=<ID=GL000225.1,length=211173,assembly=b37>
##contig=<ID=GL000192.1,length=547496,assembly=b37>
##reference=file:///cluster/project8/vyp/vincent/data/reference_genomes/fasta/human_g1k_v37.fasta
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", file = "data/sites.vcf")

write.table(x = my.vcf, file = "data/sites.vcf", append = TRUE, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)



