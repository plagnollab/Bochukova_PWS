library(ggplot2)

support <- read.table("support/Bochukova_set3.tab", header = TRUE, stringsAsFactors = FALSE)

final.frame <- data.frame()
for (i in 1:nrow(support)) {

  sample <- support$sample[ i ]
  BAM.file <- paste0("/scratch2/vyp-scratch2/Bochukova_RNASeq/processed/set3/", sample, "/", sample, "_unique.bam")
  count.file <- paste0("data/ASE/", sample, ".tab")

  counts <- read.table(count.file, sep = "\t", header = TRUE)
  counts <- subset(counts, totalCount > 20)
  counts$sample <- sample
  final.frame <- rbind.data.frame(final.frame, counts[, c("sample", "refCount", "altCount")])
  
}

### now we plot
min.freq.shown <- 0.001
my.ticks <- c(0, 0.001, 0.003, 0.01, 0.05, 0.1, 0.5, 1)

final.frame$freq <- final.frame$altCount / final.frame$refCount
final.frame$logfreq <- log10(final.frame$freq + min.freq.shown)
final.frame$depth <- final.frame$refCount + final.frame$altCount


g <- ggplot(data = final.frame, aes(x = depth, y = logfreq))
g <- g + geom_point() + facet_wrap(facets = ~ sample, nrow = 2, ncol = 4)
g <- g +  ggplot2::scale_y_continuous(breaks=log10(my.ticks + min.freq.shown), labels=my.ticks, limits = log10(range(my.ticks) + min.freq.shown))
g <- g + ggplot2::scale_x_continuous(limits = c(0, 500))
ggsave(g, file = "editing_PWS.pdf", width = 10, height = 6)
