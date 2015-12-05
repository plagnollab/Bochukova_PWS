library(pheatmap)


diff.expression.pvalues <- read.table("set3_hg38/deseq2/control_PW/deseq_Bochukova_set3_differential_expression.tab", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
diff.expression.pvalues <- diff.expression.pvalues[ !is.na(diff.expression.pvalues$external_gene_id), ]


good.genes <- diff.expression.pvalues$external_gene_id[ 1:100 ]


RPKM.values <- read.csv("/scratch2/vyp-scratch2/Bochukova_RNASeq/processed/set3_hg38/deseq2/rpkm_values.csv", stringsAsFactors = FALSE)
RPKM.values <- RPKM.values[ RPKM.values$external_gene_id %in% good.genes, ]

RPKM.mat <- RPKM.values[, grepl(pattern = "Control|PWS", names(RPKM.values)) ]
RPKM.mat <- log2(as.matrix(RPKM.mat) + 1)
dimnames(RPKM.mat)[[1]] <- RPKM.values$external_gene_id

annot.col<- data.frame(Condition = c("Control", "Control", "Control", "Control", "PWS", "PWS", "PWS", "PWS"))
row.names(annot.col) <- dimnames(RPKM.mat)[[2]]

pdf("figs/heatmap.pdf", width = 8, height = 8, onefile = FALSE)
pheatmap(mat = RPKM.mat,
         cluster_rows=TRUE,
         show_rownames=TRUE,
         cluster_cols=TRUE,
         annotation_col=annot.col,
         fontsize_row = 4)
dev.off()


