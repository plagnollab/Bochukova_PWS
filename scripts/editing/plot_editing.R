library(ggplot2)

support <- read.table("support/Bochukova_set3.tab", header = TRUE, stringsAsFactors = FALSE)

final.frame <- data.frame()
for (i in 1:nrow(support)) {

  sample <- support$sample[ i ]
  ##BAM.file <- paste0("/scratch2/vyp-scratch2/Bochukova_RNASeq/processed/set3/", sample, "/", sample, "_unique.bam")
  count.file <- paste0("data/ASE/", sample, ".tab")

  counts <- read.table(count.file, sep = "\t", header = TRUE)
  counts$ID <- paste0(counts$contig, "_", counts$position)
  
  counts <- subset(counts, totalCount > 20)
  counts$sample <- sample
  final.frame <- rbind.data.frame(final.frame, counts[, c("sample", "refCount", "altCount", "ID")])
  
}

### now we plot
min.freq.shown <- 0.001
my.ticks <- c(0, 0.001, 0.003, 0.01, 0.05, 0.1, 0.5, 1)

final.frame$freq <- final.frame$altCount / final.frame$refCount
final.frame$logfreq <- log10(final.frame$freq + min.freq.shown)
final.frame$depth <- final.frame$refCount + final.frame$altCount




if (!file.exists("figs")) dir.create("figs")

g <- ggplot(data = final.frame, aes(x = depth, y = logfreq))
g <- g + geom_point() + facet_wrap(facets = ~ sample, nrow = 2, ncol = 4)
g <- g + ggplot2::scale_y_continuous(breaks=log10(my.ticks + min.freq.shown), labels=my.ticks, limits = log10(range(my.ticks) + min.freq.shown))
g <- g + ggplot2::scale_x_continuous(limits = c(0, 500))
ggsave(g, file = "figs/editing_PWS.pdf", width = 10, height = 6)

max.freq <- tapply(final.frame$altCount / final.frame$depth, IND = final.frame$ID, FUN = max)
good.variants <- names(max.freq)[ which(max.freq > 0.05) ]


final.frame.clean <- final.frame[ final.frame$ID %in% good.variants,]
median.dispersion <- tapply(final.frame.clean$freq, FUN = median, IND = final.frame.clean$sample)
my.data <- data.frame(sample = names(median.dispersion),
                      dispersion = as.numeric(median.dispersion),
                      pheno = gsub(pattern = "[0-9]", replacement = "", names(median.dispersion)))

g <- ggplot(data = my.data, aes(x = pheno, y = dispersion)) + geom_point() + ylab(paste0("Median editing level for ", length(good.variants), "edited positions")) + xlab("Phenotype")
ggsave(g, file = "figs/editing_vs_pheno.pdf")

