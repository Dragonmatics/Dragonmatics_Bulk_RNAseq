#!/usr/bin/env Rscript

# =============================================================================
# SLURM Configuration
# =============================================================================
#SBATCH --job-name=differential_expression
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=256gb
#SBATCH --time=24:00:00
#SBATCH --output=./SLURM_OUT/%x_%A.out
#SBATCH --error=./SLURM_OUT/%x_%A.err

# =============================================================================
# Script Description
# =============================================================================
# This script performs differential gene expression analysis using DESeq2.
# It reads the gene counts data produced by the nf-core/rnaseq pipeline
# and performs pairwise comparisons between experimental conditions.

# =============================================================================
# Load Libraries
# =============================================================================
message("\n[1/13] Loading libraries........................line_23...")
library(DESeq2)
library(edgeR)
library(EnhancedVolcano)
library(ggrepel)
library(gridExtra)
library(Matrix)
library(pdftools)
library(pheatmap)
library(progress)
library(RColorBrewer)
library(tidyverse)  

# =============================================================================
# Data Import and Preprocessing
# =============================================================================
message("\n[2/13] Starting data import and preprocessing.....line_40.......")

# Read in the gene counts data
message("  Reading gene counts data...")
data <- read.table(
  "results/star_salmon/salmon.merged.gene_counts_length_scaled.tsv",
  header = TRUE,
  row.names = 1
)
data <- data[, sort(colnames(data))]
data_round <- round(data[2:9], digits = 0)
message("  Data import completed successfully")

#  Let's check the total number of reads (sequencing depth) for each sample:
# Calculate total reads per sample
total_reads <- colSums(data_round[1:8])

# Create a data frame for display
reads_df <- data.frame(
  Sample = names(total_reads),
  Total_Reads = total_reads,
  Total_Reads_Millions = round(total_reads / 1e6, 2),
  row.names = NULL
)
reads_df


# =============================================================================
# DESeq2 Object Creation and Initial Analysis
# =============================================================================
message("\n[3/13] Creating DESeq2 object and performing initial analysis.line_7.")

# Create the full DESeq2 object
message("  Setting up experimental conditions...")
condition <- c(
  rep("LNCAP_Hypoxia", 2),
  rep("LNCAP_Normoxia", 2),
  rep("PC3_Hypoxia", 2),
  rep("PC3_Normoxia", 2)
)
my_colData <- as.data.frame(condition)
my_colData$condition <- as.factor(my_colData$condition)
my_colData$rownames <- colnames(data_round)[1:8]
my_colData <- my_colData[, c("rownames", "condition")]
dds <- DESeqDataSetFromMatrix(countData = data_round,
                              colData = my_colData,
                              design = ~ condition)
message("  Running DESeq2 analysis...")
dds <- DESeq(dds)
message("  DESeq2 analysis completed successfully")

# View raw counts for first few genes
dds@assays@data$counts

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Read annotation file
annotation <- read.csv(
  "reference_genomes/GRCh38.p14_annotation.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)
# Preview annotations
head(annotation)

# Convert normalized counts to data frame with Ensembl IDs
normalized_counts <- rownames_to_column(as.data.frame(normalized_counts), var = "ensembl_id")

# Join annotations with expression data
annotated_data <- right_join(annotation,
                             normalized_counts,
                             by = c("Gene.stable.ID" = "ensembl_id"))

# Preview merged data
head(annotated_data)

# Save annotated results
write.csv(annotated_data, file = "results/gene_annotated_normalized_counts.csv", row.names = FALSE)

# Apply variance stabilizing transformation
vsd <- vst(dds, blind = TRUE)

# Confirm transformation
nrow(vsd)

# =============================================================================
# Sample Distance Analysis
# =============================================================================
message("\n[4/13] Performing sample distance analysis......line_129.........")

#' @title Plot Sample Distances
#' @description This function calculates and plots the Euclidean distance between samples.
#' @param vsd.obj A DESeqTransform object from the `vst` function.
#' @return A heatmap plot object from `pheatmap`.
plotDists <- function(vsd.obj) {
  # Calculate sample distances
  sampleDists <- dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix(sampleDists)
  
  # Add condition labels
  rownames(sampleDistMatrix) <- paste(vsd.obj$condition)
  
  # Create color palette
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
  
  # Generate heatmap
  pheatmap::pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors,
    main = "Sample Distance Matrix",
    fontsize = 15,
    border_color = NA
  )
}

# Generate and save plots
plot1 <- plotDists(vsd)
plot1

ggsave(
  filename = "figures/Euclidean_distance_between_the_expression_values.pdf",
  plot = plot1,
  width = 16,
  height = 12,
  dpi = 300
)


# =============================================================================
# Variable Gene Analysis
# =============================================================================
message("\n[5/13] Starting variable gene analysis...........line_174........")

#' @title Plot Variable Gene Heatmap
#' @description This function identifies the most variable genes and plots them as a heatmap.
#' @param vsd.obj A DESeqTransform object.
#' @param num_genes The number of top variable genes to plot.
#' @param annotation A data frame containing gene annotations.
#' @param title The title for the heatmap.
#' @return A heatmap plot object from `pheatmap`.
variable_gene_heatmap <- function(vsd.obj,
                                  num_genes = 500,
                                  annotation,
                                  title = "") {
  message(sprintf("  Analyzing top %d variable genes...", num_genes))
  # Set up color scheme
  brewer_palette <- "RdBu"
  ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]
  
  # Process expression data
  stabilized_counts <- assay(vsd.obj)
  row_variances <- rowVars(stabilized_counts)
  
  # Select top variable genes
  top_variable_genes <- stabilized_counts[order(row_variances, decreasing =
                                                  TRUE)[1:num_genes], ]
  
  # Center expression values
  top_variable_genes <- top_variable_genes -
    rowMeans(top_variable_genes, na.rm = TRUE)
  
  # Map to gene symbols
  gene_names <- annotation$Gene.name[match(rownames(top_variable_genes), annotation$Gene.stable.ID)]
  rownames(top_variable_genes) <- gene_names
  
  # Prepare sample annotations
  coldata <- as.data.frame(vsd.obj@colData)
  coldata$sizeFactor <- NULL
  
  # Generate heatmap
  pheatmap::pheatmap(
    top_variable_genes,
    color = mr,
    annotation_col = coldata,
    fontsize_col = 8,
    fontsize_row = 250 / num_genes,
    border_color = NA,
    main = title,
    show_rownames = TRUE,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    clustering_method = "complete",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean"
  )
}

# Generate and save heatmap
plot2 <- variable_gene_heatmap(vsd,
                               num_genes = 40,
                               annotation = annotation,
                               title = "Top 40 Most Variable Genes")
plot2

# Save plots
ggsave(
  filename = "figures/variable-gene-heatmap-all-samples.pdf",
  plot = plot2,
  width = 12,
  height = 9,
  dpi = 300
)

# =============================================================================
# Principal Component Analysis
# =============================================================================
message("\n[6/13] Performing Principal Component Analysis........line_250......")

#' @title Plot PCA
#' @description This function performs PCA and generates a scatter plot.
#' @param vsd.obj A DESeqTransform object.
#' @return A ggplot object.
plot_PCA <- function(vsd.obj) {
  # Generate PCA data
  pcaData <- plotPCA(vsd.obj,
                     intgroup = c("condition"),
                     returnData = TRUE)
  
  # Calculate variance percentages
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  # Create enhanced PCA plot
  ggplot(pcaData, aes(PC1, PC2, color = condition)) +
    # Add points with larger size
    geom_point(size = 8, alpha = 0.8) +
    
    # Add labels with repulsion
    ggrepel::geom_text_repel(
      aes(label = name),
      color = "black",
      size = 3,
      box.padding = 0.5,
      point.padding = 0.5
    ) +
    
    # Customize theme
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5
      )
    ) +
    
    # Add labels and title
    labs(
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      title = "Principal Component Analysis",
      color = "Condition"
    )
}

# Generate and save PCA plot
plot3 <- plot_PCA(vsd)
plot3

ggsave(
  filename = "figures/PCA-plot.pdf",
  plot = plot3,
  width = 16,
  height = 12,
  dpi = 300
)


# =============================================================================
# Cell-Specific Analysis
# =============================================================================
message("\n[7/13] Starting cell-specific analysis.....line_321............")

#' @title Generate DESeq2 Object
#' @description Creates a DESeq2 object for a subset of data.
#' @param my_data The full count data matrix.
#' @param groups A character vector with two group names to compare.
#' @return A DESeqDataSet object.
generate_DESeq_object <- function (my_data, groups) {
  data_subset1 <- my_data[, grep(str_c("^", groups[1]), colnames(my_data))]
  data_subset2 <- my_data[, grep(str_c("^", groups[2]), colnames(my_data))]
  my_countData <- cbind(data_subset1, data_subset2)
  condition <- c(rep(groups[1], ncol(data_subset1)), rep(groups[2], ncol(data_subset2)))
  my_colData <- as.data.frame(condition)
  rownames(my_colData) <- colnames(my_countData)
  dds <- DESeqDataSetFromMatrix(countData = my_countData,
                                colData = my_colData,
                                design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)
  return(dds)
}

# Create subsetted DESeq2 objects for LNCaP and PC3 cells
lncap <- generate_DESeq_object(data_round, c("LNCAP_Hypoxia", "LNCAP_Normoxia"))
pc3 <- generate_DESeq_object(data_round, c("PC3_Hypoxia", "PC3_Normoxia"))

# Perform variance stabilizing transformation
lncap_vsd <- vst(lncap, blind = TRUE)
pc3_vsd <- vst(pc3, blind = TRUE)

# Generate LNCaP variable genes heatmap
plot4 <- variable_gene_heatmap(
  lncap_vsd,
  num_genes = 30,
  annotation = annotation,
  title = "LNCaP Hypoxia Response Genes"
)
# Save cell-specific heatmaps
ggsave(
  filename = "figures/LNCaP_variable_genes.pdf",
  plot = plot4,
  width = 12,
  height = 8
)

# Generate PC3 variable genes heatmap

plot5 <- variable_gene_heatmap(pc3_vsd, 30, annotation = annotation, title = "PC3 variable genes")
pdf(file = "figures/PC3_LNCaP_variable_genes.pdf",
    width = 16,
    height = 12)
plot6 <- gridExtra::grid.arrange(plot4[[4]], plot5[[4]], nrow = 1)
dev.off()


# Extract differentially expressed genes between hypoxia and normoxia conditions using DESeq2's `results()` function:

results(lncap,
        contrast = c("condition", "LNCAP_Hypoxia", "LNCAP_Normoxia"))

head(as.data.frame(results(
  lncap,
  contrast = c("condition", "LNCAP_Hypoxia", "LNCAP_Normoxia")
)))


# =============================================================================
# Differential Expression Analysis
# =============================================================================
message("\n[8/13] Starting cell-specific analysis.....line_389............")

#' @title Generate Differential Expression Results
#' @description This function performs DE analysis and saves the results to CSV files.
#' @param dds A DESeqDataSet object.
#' @param comparisons A character vector with two group names to compare.
#' @param padjcutoff Adjusted p-value cutoff for significance.
#' @param log2cutoff Log2 fold change cutoff for significance.
#' @param cpmcutoff CPM cutoff for filtering.
#' @return A tibble summarizing the number of significant genes.
generate_DE_results <- function (dds,
                                 comparisons,
                                 padjcutoff = 0.001,
                                 log2cutoff = 0.5,
                                 cpmcutoff = 2) {
  message(
    sprintf(
      "\nAnalyzing differential expression for %s vs %s:",
      comparisons[1],
      comparisons[2]
    )
  )
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent eta: :eta",
    total = 7,
    # Number of major steps in this function
    width = 60
  )
  raw_counts <- counts(dds, normalized = F)
  pb$tick()
  cpms <- enframe(rowMeans(edgeR::cpm(raw_counts)))
  colnames(cpms) <- c("ensembl_id", "avg_cpm")
  pb$tick()
  res <- results(dds, contrast = c("condition", comparisons[1], comparisons[2]))[, -c(3, 4)]
  res <- as_tibble(res, rownames = "ensembl_id")
  my_annotation <- read.csv(
    "reference_genomes/GRCh38.p14_annotation.csv",
    header = TRUE,
    stringsAsFactors = F
  )
  res <- left_join(res, my_annotation, by = c("ensembl_id" = "Gene.stable.ID"))
  res <- left_join(res, cpms, by = c("ensembl_id" = "ensembl_id"))
  normalized_counts <- round(counts(dds, normalized = TRUE), 3)
  pattern <- str_c(comparisons[1], "|", comparisons[2])
  combined_data <- as_tibble(cbind(res, normalized_counts[, grep(pattern, colnames(normalized_counts))]))
  combined_data <- combined_data[order(combined_data$log2FoldChange, decreasing = TRUE), ]
  res_prot <- res[which(res$Gene.type == "protein_coding"), ]
  res_prot_ranked <- res_prot[order(res_prot$log2FoldChange, decreasing = TRUE), c("Gene.name", "log2FoldChange")]
  res_prot_ranked <- na.omit(res_prot_ranked)
  res_prot_ranked <- res_prot_ranked[res_prot_ranked$Gene.name != "", ]
  res_prot_ranked$Gene.name <- str_to_upper(res_prot_ranked$Gene.name)
  de_genes_padj <- res[which(res$padj < padjcutoff), ]
  de_genes_log2f <- res[which(abs(res$log2FoldChange) > log2cutoff &
                                res$padj < padjcutoff), ]
  de_genes_cpm <- res[which(res$avg_cpm > cpmcutoff &
                              res$padj < padjcutoff), ]
  write.csv (
    de_genes_padj,
    file = paste0(
      "results/",
      comparisons[1],
      "_vs_",
      comparisons[2],
      "_padj_cutoff.csv"
    ),
    row.names = FALSE
  )
  write.csv (
    de_genes_log2f,
    file = paste0(
      "results/",
      comparisons[1],
      "_vs_",
      comparisons[2],
      "_log2f_cutoff.csv"
    ),
    row.names = FALSE
  )
  write.csv (
    de_genes_cpm,
    file = paste0(
      "results/",
      comparisons[1],
      "_vs_",
      comparisons[2],
      "_cpm_cutoff.csv"
    ),
    row.names = FALSE
  )
  write.csv (
    combined_data,
    file = paste0(
      "results/",
      comparisons[1],
      "_vs_",
      comparisons[2],
      "_allgenes.csv"
    ),
    row.names = FALSE
  )
  write.table (
    res_prot_ranked,
    file = paste0("results/", comparisons[1], "_vs_", comparisons[2], "_rank.rnk"),
    sep = "\t",
    row.names = F,
    quote = F
  )
  
  writeLines(
    paste0(
      "For the comparison: ",
      comparisons[1],
      "_vs_",
      comparisons[2],
      ", out of ",
      nrow(combined_data),
      " genes, there were: \n",
      nrow(de_genes_padj),
      " genes below padj ",
      padjcutoff,
      "\n",
      nrow(de_genes_log2f),
      " genes below padj ",
      padjcutoff,
      " and above a log2FoldChange of ",
      log2cutoff,
      "\n",
      nrow(de_genes_cpm),
      " genes below padj ",
      padjcutoff,
      " and above an avg cpm of ",
      cpmcutoff,
      "\n",
      "Gene lists ordered by log2fchange with the cutoffs above have been generated."
    )
  )
  gene_count <- tibble (
    cutoff_parameter = c("padj", "log2fc", "avg_cpm"),
    cutoff_value = c(padjcutoff, log2cutoff, cpmcutoff),
    signif_genes = c(
      nrow(de_genes_padj),
      nrow(de_genes_log2f),
      nrow(de_genes_cpm)
    )
  )
  invisible(gene_count)
}

# Generate and save the DE results
lncap_output <- generate_DE_results (lncap, c("LNCAP_Hypoxia", "LNCAP_Normoxia"))
pc3_output <- generate_DE_results(pc3, c("PC3_Hypoxia", "PC3_Normoxia"))

res <- read.csv("results/LNCAP_Hypoxia_vs_LNCAP_Normoxia_allgenes.csv",
                header = TRUE)
head(res)

# =============================================================================
# Results Visualization
# =============================================================================
message("\n[9/13] Generating visualization plots........line_548............")

# Function for plotting individual gene counts
d <- plotCounts(dds,
                gene = "ENSG00000146678",
                intgroup = "condition",
                returnData = TRUE)

plot7 <- ggplot(d, aes(x = condition, y = count)) +
  geom_point(size = 10,
             position = position_jitter(w = 0.4, h = 0),
             aes (color = condition)) +   scale_y_log10(breaks = c(25, 100, 400))
plot7

ggsave(
  filename = "figures/ENSG00000146678.pdf",
  plot = plot7,
  width = 16,
  height = 12
)


#' @title Plot Gene Counts
#' @description This function plots the normalized counts for a specific gene.
#' @param dds A DESeqDataSet object.
#' @param gene The Ensembl ID or gene name to plot.
#' @param normalization The normalization method to use ("cpm" or "DESeq2").
#' @return A ggplot object.
plot_counts <- function (dds, gene, normalization = "DESeq2") {
  # read in the annotation file
  annotation <- read.csv(
    "reference_genomes/GRCh38.p14_annotation.csv",
    header = TRUE,
    stringsAsFactors = F
  )
  # obtain normalized data
  if (normalization == "cpm") {
    normalized_data <- cpm(counts(dds, normalized = F)) # normalize the raw data by counts per million
  } else if (normalization == "DESeq2")
    normalized_data <- counts(dds, normalized = TRUE) # use DESeq2 normalized counts
  # get sample groups from colData
  condition <- dds@colData$condition
  # get the gene name from the ensembl id
  if (is.numeric(gene)) {
    # check if an index is supplied or if ensembl_id is supplied
    if (gene %% 1 == 0)
      ensembl_id <- rownames(normalized_data)[gene]
    else
      stop("Invalid index supplied.")
  } else if (gene %in% annotation$Gene.name) {
    # check if a gene name is supplied
    ensembl_id <- annotation$Gene.stable.ID[which(annotation$Gene.name == gene)]
  } else if (gene %in% annotation$Gene.stable.ID) {
    ensembl_id <- gene
  } else {
    stop("Gene not found. Check spelling.")
  }
  expression <- normalized_data[ensembl_id, ]
  gene_name <- annotation$Gene.name[which(annotation$Gene.stable.ID == ensembl_id)]
  # construct a tibble with the grouping and expression
  gene_tib <- tibble(condition = condition, expression = expression)
  ggplot(gene_tib, aes(x = condition, y = expression)) +
    geom_boxplot(aes(fill = condition), outlier.size = NULL) +
    geom_point(aes(color = condition)) +
    labs (
      title = paste0("Expression of ", gene_name, " - ", ensembl_id),
      x = "group",
      y = paste0("Normalized expression (", normalization , ")")
    ) +
    theme(axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11))
}

plot8 <- plot_counts(dds, "IGFBP1")
plot8

ggsave(
  filename = "figures/IGFBP1_plot_counts.pdf",
  plot = plot8,
  width = 16,
  height = 12
)


# =============================================================================
# Differential Gene Heatmap
# =============================================================================
message("\n[10/13] Differential Gene Heatmap........line_635...............")
res <- read.csv ("results/LNCAP_Hypoxia_vs_LNCAP_Normoxia_allgenes.csv",
                 header = TRUE)

#' @title Plot DE Gene Heatmap
#' @description This function plots a heatmap of differentially expressed genes.
#' @param res A data frame with DE results.
#' @param padj_cutoff Adjusted p-value cutoff for significance.
#' @param ngenes The number of genes to include in the heatmap.
#' @return A pheatmap object.
DE_gene_heatmap <- function(res,
                            padj_cutoff = 0.0001,
                            ngenes = 20) {
  # generate the color palette
  brewer_palette <- "RdBu"
  ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]
  # obtain the significant genes and order by log2FoldChange
  significant_genes <- res %>% filter(padj < padj_cutoff) %>% arrange (desc(log2FoldChange)) %>% head (ngenes)
  heatmap_values <- as.matrix(significant_genes[, -c(1:8)])
  rownames(heatmap_values) <- significant_genes$Gene.name
  # plot the heatmap using pheatmap
  pheatmap::pheatmap(
    heatmap_values,
    color = mr,
    scale = "row",
    fontsize_col = 12,
    fontsize_row = 200 / ngenes,
    fontsize = 5,
    border_color = NA
  )
}

plot9 <- DE_gene_heatmap(res, 0.001, 30)

ggsave(filename = "figures/DE_gene_heatmap.pdf",
       plot = plot9,
       width = 16,
       height = 12)


# =============================================================================
# Plot Volcano
# =============================================================================
message("\n[11/13] Plot Volcano........................line_635...............")


#' @title Plot Volcano
#' @description This function generates a volcano plot from DE results.
#' @param res A data frame with DE results.
#' @param padj_cutoff Adjusted p-value cutoff for significance.
#' @param fc_cutoff Fold change cutoff for significance.
#' @param nlabel The number of genes to label.
#' @param label.by The metric to use for labeling top genes ("padj" or "log2FoldChange").
#' @return A ggplot object.
plot_volcano <- function(res,
                         padj_cutoff = 0.05,
                         fc_cutoff = 1,
                         nlabel = 10,
                         label.by = "padj") {
  # Create significance categories
  res <- res %>%
    mutate(
      expression = case_when(
        log2FoldChange >= fc_cutoff & padj < padj_cutoff ~ "Upregulated",
        log2FoldChange <= -fc_cutoff &
          padj < padj_cutoff ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    ) %>%
    filter(!is.na(padj))  # Remove NA values
  
  # Prepare genes for labeling
  if (label.by == "padj") {
    up_genes <- res %>%
      filter(expression == "Upregulated") %>%
      arrange(padj) %>%
      head(nlabel)
    down_genes <- res %>%
      filter(expression == "Downregulated") %>%
      arrange(padj) %>%
      head(nlabel)
  } else if (label.by == "log2FoldChange") {
    up_genes <- res %>%
      filter(expression == "Upregulated") %>%
      arrange(desc(abs(log2FoldChange))) %>%
      head(nlabel)
    down_genes <- res %>%
      filter(expression == "Downregulated") %>%
      arrange(desc(abs(log2FoldChange))) %>%
      head(nlabel)
  }
  
  # Generate statistics for title
  stats <- list(
    total = nrow(res),
    up = sum(res$expression == "Upregulated"),
    down = sum(res$expression == "Downregulated")
  )
  
  # Create enhanced volcano plot
  ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
    # Add points
    geom_point(aes(color = expression), alpha = 0.7, size = 1) +
    # Color scheme
    scale_color_manual(
      values = c(
        "Upregulated" = "#fc8d59",
        "Downregulated" = "#91bfdb",
        "Not Significant" = "grey80"
      )
    ) +
    # Add gene labels with repulsion
    ggrepel::geom_text_repel(
      data = up_genes,
      aes(label = Gene.name),
      size = 3,
      color = "#d73027",
      max.overlaps = Inf,
      box.padding = 0.5
    ) +
    ggrepel::geom_text_repel(
      data = down_genes,
      aes(label = Gene.name),
      size = 3,
      color = "#4575b4",
      max.overlaps = Inf,
      box.padding = 0.5
    ) +
    # Add threshold lines
    geom_vline(
      xintercept = c(-fc_cutoff, fc_cutoff),
      linetype = "dashed",
      color = "grey50",
      alpha = 0.5
    ) +
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed",
      color = "grey50",
      alpha = 0.5
    ) +
    # Customize theme
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 12),
      plot.subtitle = element_text(size = 10)
    ) +
    # Add labels
    labs(
      title = "Differential Expression Analysis",
      subtitle = sprintf(
        "Total: %d genes | Up: %d | Down: %d (padj < %.2g, |log2FC| > %.1f)",
        stats$total,
        stats$up,
        stats$down,
        padj_cutoff,
        fc_cutoff
      ),
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)",
      color = "Expression Change"
    )
}

# Generate volcano plots with different thresholds and labeling strategies
#  Default stringent analysis
plot10 <- plot_volcano(
  res,
  padj_cutoff = 0.001,
  # Stringent significance threshold
  fc_cutoff = 1.5,
  # Require at least 2.8-fold change
  nlabel = 10,
  # Label top 10 genes
  label.by = "padj"
)    # Label most significant genes

plot10

ggsave(
  filename = "figures/Volcano_plot_stringent.pdf",
  plot = plot10,
  width = 12,
  height = 10
)


#  Alternative analysis focusing on fold change
plot11 <- plot_volcano(
  res,
  padj_cutoff = 0.05,
  # More permissive significance
  fc_cutoff = 2,
  # Focus on genes with larger changes
  nlabel = 15,
  # Label more genes
  label.by = "log2FoldChange"
)  # Label genes with largest changes

plot11

ggsave(
  filename = "figures/Volcano_plot.pdf",
  plot = plot11,
  width = 12,
  height = 10
)


# =============================================================================
# Enhanced volcano plots
# =============================================================================
message("\n[12/13] Plot Volcano........................line_635...............")

# Create a column for gene names from the rownames
res$gene_names <- res$Gene.name

# Sort the results by padj and then by log2FoldChange
# This is a good way to find the most significant and most-changed genes
sorted_res <- res[order(res$padj, abs(res$log2FoldChange), decreasing = c(FALSE, TRUE)), ]

# Filter for significantly upregulated genes
top_up_genes <- subset(sorted_res, log2FoldChange > 2 & padj < 0.05)

# Filter for significantly downregulated genes
top_down_genes <- subset(sorted_res, log2FoldChange < -2 & padj < 0.05)

# Select the top 10 gene names from each list
select_lab_list <- c(head(top_up_genes, 10)$gene_names, head(top_down_genes, 10)$gene_names)

# Create the enhanced volcano plot
plot12 <- EnhancedVolcano(res,
                       lab = res$gene_name, # Use the column with gene names
                       x = 'log2FoldChange',
                       y = 'padj',
                       ylab = bquote(~-Log[10] ~ italic(FDR)),
                       pCutoff = 0.05,
                       FCcutoff = 2,
                       selectLab = select_lab_list, # Add this line
                       legendLabels=c('Not DE &\nabsolute FC < 2',
                                      'Not DE &\nabsolute FC > 2',
                                      'FDR < 0.05 &\nabsolute FC < 2',
                                      'FDR < 0.05 &\nabsolute FC > 2'
                       ),
                       legendPosition = 'right')

plot12

ggsave(
  filename = paste0("figures/Volcanoplot-enhanced.pdf"),
  plot = plot12,
  width = 7,
  height = 7
)

# Save results with detailed categorization
res_processed <- res %>%
  mutate(
    expression_category = case_when(
      log2FoldChange >= 2 & padj <= 0.001 ~ "Strongly Upregulated",
      log2FoldChange >= 1 & padj <= 0.05 ~ "Upregulated",
      log2FoldChange <= -2 &
        padj <= 0.001 ~ "Strongly Downregulated",
      log2FoldChange <= -1 & padj <= 0.05 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

# Save results for downstream analysis
saveRDS(res_processed,
        "results/processed_differential_expression.rds")
write.csv(res_processed,
          "results/processed_differential_expression.csv",
          row.names = FALSE)



# =============================================================================
# LogFoldChange comparison plot
# =============================================================================
message("\n[13/13] LogFoldChange comparison plot.........line_918...............")

res1 <- read.csv ("results/LNCAP_Hypoxia_vs_LNCAP_Normoxia_allgenes.csv",
                  header = TRUE)
res2 <- read.csv ("results/PC3_Hypoxia_vs_PC3_Normoxia_allgenes.csv", header = TRUE)

#' @title Compare Significant Genes
#' @description This function compares significant genes between two result sets and plots the log2FoldChange values against each other.
#' @param res1 A data frame with DE results for the first comparison.
#' @param res2 A data frame with DE results for the second comparison.
#' @param padj_cutoff Adjusted p-value cutoff for significance.
#' @param ngenes The number of top up- and down-regulated genes to consider from each result set.
#' @param nlabel The number of genes to label in the plot.
#' @param samplenames A character vector with two sample names for the plot axes.
#' @param title The title of the plot.
#' @return A ggplot object.
compare_significant_genes <- function (res1,
                                       res2,
                                       padj_cutoff = 0.0001,
                                       ngenes = 250,
                                       nlabel = 10,
                                       samplenames = c("comparison1", "comparison2"),
                                       title = "") {
  # get list of most upregulated or downregulated genes for each results table
  genes1 <- rbind(head(res1[which(res1$padj < padj_cutoff), ], ngenes), tail(res1[which(res1$padj < padj_cutoff), ], ngenes))
  genes2 <- rbind(head(res2[which(res2$padj < padj_cutoff), ], ngenes), tail(res2[which(res2$padj < padj_cutoff), ], ngenes))
  
  # combine the data from both tables
  de_union <- union(genes1$ensembl_id, genes2$ensembl_id)
  res1_union <- res1[match(de_union, res1$ensembl_id), ][c("ensembl_id", "log2FoldChange", "Gene.name")]
  res2_union <- res2[match(de_union, res2$ensembl_id), ][c("ensembl_id", "log2FoldChange", "Gene.name")]
  combined <- left_join(res1_union, res2_union, by = "ensembl_id", suffix = samplenames)
  
  # identify overlap between genes in both tables
  combined$de_condition <- 1 # makes a placeholder column
  combined$de_condition[which(combined$ensembl_id %in% intersect(genes1$ensembl_id, genes2$ensembl_id))] <- "Significant in Both"
  combined$de_condition[which(combined$ensembl_id %in% setdiff(genes1$ensembl_id, genes2$ensembl_id))] <- paste0("Significant in ", samplenames[1])
  combined$de_condition[which(combined$ensembl_id %in% setdiff(genes2$ensembl_id, genes1$ensembl_id))] <- paste0("Significant in ", samplenames[2])
  combined[is.na(combined)] <- 0
  
  # find the top most genes within each condition to label on the graph
  label1 <- rbind(head(combined[which(combined$de_condition == paste0("Significant in ", samplenames[1])), ], nlabel), tail(combined[which(combined$de_condition ==
                                                                                                                                             paste0("Significant in ", samplenames[1])), ], nlabel))
  label2 <- rbind(head(combined[which(combined$de_condition == paste0("Significant in ", samplenames[2])), ], nlabel), tail(combined[which(combined$de_condition ==
                                                                                                                                             paste0("Significant in ", samplenames[2])), ], nlabel))
  label3 <- rbind(head(combined[which(combined$de_condition == "Significant in Both"), ], nlabel), tail(combined[which(combined$de_condition ==
                                                                                                                         "Significant in Both"), ], nlabel))
  combined_labels <- rbind(label1, label2, label3)
  
  # plot the genes based on log2FoldChange, color coded by significance
  ggplot(combined, aes_string(
    x = paste0("log2FoldChange", samplenames[1]),
    y = paste0("log2FoldChange", samplenames[2])
  )) +
    geom_point(aes(color = de_condition), size = 1.7) +
    scale_color_manual(values = c("#00BA38", "#619CFF", "#F8766D", "red")) +
    ggrepel::geom_text_repel(
      data = combined_labels,
      aes_string(
        label = paste0("Gene.name", samplenames[1]),
        color = "de_condition"
      ),
      show.legend = F,
      size = 3
    ) +
    geom_vline(xintercept = c(0, 0),
               size = 0.3,
               linetype = 2) +
    geom_hline(yintercept = c(0, 0),
               size = 0.3,
               linetype = 2) +
    labs(
      title = title,
      x = paste0("log2FoldChange in ", samplenames[1]),
      y = paste0("log2FoldChange in ", samplenames[2])
    ) +
    theme_minimal() +
    theme(legend.title = element_blank())
}

plot13 <- compare_significant_genes(res1,
                                 res2,
                                 samplenames = c("LNCaP", "PC3"),
                                 title = "Hypoxia-induced gene expression differences in LNCaP vs PC3 cells")

plot13

ggsave(
  filename = "figures/compare_significant_genes.pdf",
  plot = plot13,
  width = 16,
  height = 12
)


# Save the DESeq2 objects for later use
saveRDS(dds, file = "results/dds.rds")
saveRDS(lncap, file = "results/lncap_dds.rds")
saveRDS(pc3, file = "results/pc3_dds.rds")

print("Differential expression analysis complete.")
