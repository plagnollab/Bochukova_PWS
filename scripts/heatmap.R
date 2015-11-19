library(pheatmap)
library(DESeq2)


load("set3_hg38/deseq2/control_PW/deseq2_object.RData")


annotations <- read.table("/scratch2/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/biomart_annotations_human.tab", header = TRUE, stringsAsFactors = FALSE)




my.results <- results(CDS)

pdf("figs/MAplot.pdf")
plotMA(my.results, main="MA-plot", ylim=c(-4,4))
dev.off()



select <-   order(my.results$pvalue)[1:100]
nt <- rlog(CDS) # defaults to log2(x+1)


pdf("figs/PCA.pdf")
plotPCA(nt, intgroup=c("condition"))
dev.off()


log2.norm.counts <- assay(nt)[select,]

clean.names <- annotations$external_gene_name[ match(rownames(log2.norm.counts), annotations$EnsemblID)  ]
good.names <- !is.na(clean.names)
log2.norm.counts <- log2.norm.counts[ good.names, ]
rownames(log2.norm.counts) <- clean.names [ good.names ]

df <- as.data.frame(colData(CDS)[,c("condition")])
names(df) <- "Condition"


pdf("figs/heatmap.pdf", width = 8, height = 8, onefile = FALSE)
pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, fontsize_row = 4)
dev.off()
