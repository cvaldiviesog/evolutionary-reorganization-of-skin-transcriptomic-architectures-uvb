# Supplementary Text S1
# R script for differential gene expression analysis
#
# The same analytical workflow was applied to O. ascotanensis, wild-type
# zebrafish, and nacre zebrafish. Only the input count matrix, sample
# identifiers, and experimental group assignments differed among datasets.

library(edgeR)
library(ggplot2)

# -------------------------------------------------------------------------
# Differential expression analysis function
# -------------------------------------------------------------------------
#
# Arguments:
#   count_file      tab-delimited count matrix
#   group_file      tab-delimited sample metadata
#   sample_columns  columns in the count matrix corresponding to samples
#   group_vector    factor defining the experimental groups
#   gene_column     column containing gene identifiers
#   output_prefix   prefix for output files
#
# Differentially expressed genes are selected using:
#   |logFC| > 1
#   P-value < 0.001
#   FDR < 0.05
#
# The same workflow is applied to all three datasets.

run_DE_analysis <- function(count_file,
                            group_file,
                            sample_columns,
                            group_vector,
                            gene_column = 1,
                            output_prefix) {

  counts <- read.delim(
    file = count_file,
    sep = "\t",
    header = TRUE,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  counts <- na.omit(counts)

  gene_ids <- counts[[gene_column]]

  count_matrix <- counts[, sample_columns, drop = FALSE]
  rownames(count_matrix) <- make.names(gene_ids, unique = TRUE)

  group_vector <- factor(group_vector)

  dge <- DGEList(
    count_matrix,
    group = group_vector,
    remove.zeros = TRUE,
    genes = data.frame(GeneID = gene_ids)
  )

  # Exploratory assessment of count distributions
  pseudo_counts <- log2(dge$counts + 1)

  pdf(paste0(output_prefix, ".counts_distribution.pdf"))
  boxplot(
    pseudo_counts,
    las = 3,
    cex.names = 1,
    main = "Log2-transformed raw counts"
  )
  dev.off()

  # Dispersion estimation
  dge <- estimateCommonDisp(dge)
  dge <- estimateTagwiseDisp(dge)

  pdf(paste0(output_prefix, ".BCV.pdf"))
  plotBCV(dge)
  dev.off()

  # Multidimensional scaling
  pdf(paste0(output_prefix, ".MDS.pdf"))
  plotMDS.DGEList(dge)
  dev.off()

  # TMM normalization
  dge <- calcNormFactors(dge, method = "TMM")

  normalized_counts <- cpm(dge)
  log_normalized_counts <- cpm(
    dge,
    log = TRUE,
    prior.count = 1
  )

  pdf(paste0(output_prefix, ".normalized_counts.pdf"))
  boxplot(
    log_normalized_counts,
    las = 3,
    cex.names = 1,
    main = "Log2-CPM normalized counts"
  )
  dev.off()

  # Differential expression testing
  exact_test <- exactTest(dge)

  results <- topTags(
    exact_test,
    n = nrow(exact_test$table)
  )$table

  # Differentially expressed genes
  is_DE <- results$FDR < 0.05 &
           abs(results$logFC) > 1 &
           results$PValue < 0.001

  DE_genes <- results[is_DE, , drop = FALSE]

  write.table(
    DE_genes,
    file = paste0(output_prefix, ".DE_genes.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE
  )

  # Regulation status for visualization and supplementary tables
  results$DE <- "NS"
  results$DE[
    results$logFC > 1 &
    results$PValue < 0.001
  ] <- "UP"

  results$DE[
    results$logFC < -1 &
    results$PValue < 0.001
  ] <- "DOWN"

  write.table(
    results,
    file = paste0(output_prefix, ".all_genes.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = TRUE
  )

  # Volcano plot
  pdf(paste0(output_prefix, ".volcano.pdf"))

  ggplot(
    results,
    aes(
      x = logFC,
      y = -log10(PValue),
      colour = DE
    )
  ) +
    geom_point() +
    theme_classic() +
    labs(
      x = "logFC",
      y = "-log10(p-value)",
      colour = "DE status"
    )

  dev.off()

  return(
    list(
      DGEList = dge,
      results = results,
      DE_genes = DE_genes,
      normalized_counts = normalized_counts
    )
  )
}

