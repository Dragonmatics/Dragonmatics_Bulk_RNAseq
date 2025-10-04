#!/usr/bin/env Rscript

# Set the start time for the entire script
Sys.setenv(START_TIME = format(Sys.time()))

# =============================================================================
# SLURM Configuration
# =============================================================================
#SBATCH --job-name=functional_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=256gb
#SBATCH --time=22:00:00
#SBATCH --output=./SLURM_OUT/%x_%A_%a.out 
#SBATCH --error=./SLURM_OUT/%x_%A_%a.err 

# =============================================================================
# Script Description
# =============================================================================
# This script performs functional analysis, including GSEA, GO, and KEGG pathway analysis.
# It uses the results from the differential expression analysis to identify enriched
# biological pathways and functions.

# =============================================================================
# Load Required Libraries
# =============================================================================
# Function to print timestamped messages
print_status <- function(msg) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(sprintf("[%s] %s", timestamp, msg))
}

print_status("Starting script execution...")
print_status("Loading required libraries...")

# Create a function to track package loading
load_package <- function(package_name) {
    print_status(sprintf("Loading %s...", package_name))
    library(package_name, character.only = TRUE)
}

# Essential packages for analysis
load_package( "biomaRt") # For gene annotations
load_package( "circlize") # For Circos plots
load_package( "clusterProfiler") # For enrichment analysis
load_package( "DESeq2") # For differential expression results
load_package( "enrichplot" ) # For enrichment visualization
load_package( "fgsea" ) # For GSEA analysis
load_package( "ggplot2" ) # For plotting
load_package( "org.Hs.eg.db" ) # For gene ID conversion
# load_package( "pathview") # For pathway visualization
load_package( "tidyverse") # For data manipulation

# =============================================================================
# Setup Output Directories
# =============================================================================
print_status("\1 Creating output directories........line_57................\1")
dir.create("pathview_figures", showWarnings = FALSE)
dir.create("SBGNview_figures", showWarnings = FALSE)

# =============================================================================
# GSEA Analysis
# =============================================================================
print_status("\2Starting GSEA Analysis................line_64...........\2")
# Load hallmark pathways
print_status("\4Loading hallmark pathways................line_66...........\4")
hallmark_pathway <- gmtPathways("reference_genomes/h.all.v2024.1.Hs.symbols.gmt.txt")
head(names(hallmark_pathway))

print_status("\3Preparing ranked gene list................line_70...........\3")
prepare_ranked_list <- function(ranked_list) { 
  if( sum(duplicated(ranked_list$Gene.name)) > 0) {
    ranked_list <- aggregate(.~Gene.name, FUN = mean, data = ranked_list)
    ranked_list <- ranked_list[order(ranked_list$log2FoldChange, decreasing = TRUE),]
  }
  ranked_list <- na.omit(ranked_list)
  ranked_list <- tibble::deframe(ranked_list)
  ranked_list
}

# load the ranked list
lncap_ranked_list <- read.table("results/LNCAP_Hypoxia_vs_LNCAP_Normoxia_rank_2.rnk", header = TRUE, stringsAsFactors = F)
lncap_ranked_list <- prepare_ranked_list(lncap_ranked_list)

# Run fgsea
print_status("\5Running GSEA analysis (this may take a while)................line_85...........\5")
start_time <- Sys.time()
fgsea_results <- fgsea(pathways = hallmark_pathway,
                  stats = lncap_ranked_list,
                  minSize = 15,
                  maxSize = 500,
                  nperm= 1000)
end_time <- Sys.time()
print_status(sprintf("GSEA analysis completed in %s", format(end_time - start_time)))

fgsea_results %>% arrange (desc(NES)) %>% dplyr::select (pathway, padj, NES) %>% head()


# Waterfall plot
waterfall_plot <- function (fsgea_results, graph_title) {
  fgsea_results %>% 
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
    ggplot( aes(reorder(short_name,NES), NES)) +
      geom_bar(stat= "identity", aes(fill = padj<0.05))+
      coord_flip()+
      labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
      theme(axis.text.y = element_text(size = 7), 
            plot.title = element_text(hjust = 1))
}
p1 <- waterfall_plot(fgsea_results, "Hallmark pathways altered by hypoxia in LNCaP cells")
p1
ggsave(filename = "figures/waterfall_plot.pdf", plot = p1, width = 16, height = 12)



# GO and KEGG analysis
# Prepare input for clusterProfiler
df = read.csv("results/LNCAP_Hypoxia_vs_LNCAP_Normoxia_allgenes.csv")
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$ensembl_id
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

# wrapper for fgsea::plotEnrichment()
plot_enrichment <- function (geneset, pathway, ranked_list) {
  plotEnrichment(geneset[[pathway]], ranked_list)+labs (title = pathway)
}
# example of positively enriched pathway (up in Hypoxia)
p2 <- plot_enrichment(hallmark_pathway, "HALLMARK_GLYCOLYSIS" , lncap_ranked_list)
p2

ggsave(filename = "figures/HALLMARK_GLYCOLYSIS.pdf", plot = p2, width = 16, height = 12)


p3 <- plot_enrichment(hallmark_pathway, "HALLMARK_OXIDATIVE_PHOSPHORYLATION" , lncap_ranked_list)
p3

ggsave(filename = "figures/HALLMARK_OXIDATIVE_PHOSPHORYLATION.pdf", plot = p3, width = 16, height = 12)



# =============================================================================
# GO Analysis: Data Preparation
# =============================================================================
print_status("\6Starting GO analysis preparation................line_157...........\6")
print_status("\7Reading DESeq2 results................line_158...........\7")
df = read.csv("results/res_2.csv")
# df = read.csv("data/drosphila_example_de.csv", header=TRUE)
# we want the log2 fold change 
original_gene_list <- df$log2FoldChange
# name the vector
names(original_gene_list) <- df$ensembl_id
# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
head(gene_list)

# =============================================================================
# Gene Set Enrichment Analysis
# =============================================================================
print_status("\1Starting Gene Set Enrichment Analysis................line_177...........\1")
start_time <- Sys.time()
gse <- clusterProfiler::gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

# GO plots
p4 <- dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
p4
ggsave(filename = "figures/gse_dotplot.pdf", plot = p4, width = 16, height = 12)




gse_2 <- pairwise_termsim(gse, showCategory = dim(gse)[1] )

geneSets2plot <- p4$data$Description
p5 <- enrichplot::emapplot(gse_2, showCategory=geneSets2plot)
p5
ggsave(filename = "figures/gse_emapplot.pdf", plot = p5, width = 16, height = 12)


### 🌳 Tree plot
p6 <- treeplot(gse_2)
p6

ggsave(filename = "figures/gse_treeplot.pdf", plot = p6, width = 16, height = 12)


### 🕸️ Category Netplot
# categorySize can be either 'pvalue' or 'geneNum'
p7 <- cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
p7

ggsave(filename = "figures/gse_cnetplot.pdf", plot = p7, width = 16, height = 12)


### 📊 UpSet Plot
p8 <- upsetplot(gse)
p8

ggsave(filename = "figures/gse_upsetplot.pdf", plot = p8, width = 16, height = 12)


### ⛰️ Ridgeplot
p9 <- ridgeplot(gse) + labs(x = "enrichment distribution")
p9

ggsave(filename = "figures/gse_redgeplot.pdf", plot = p9, width = 16, height = 12)


### 📊 GSEA Plot
p10 <- gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
p10

ggsave(filename = "figures/gse_gseaplot.pdf", plot = p10, width = 16, height = 12)


### 📚 PubMed trend of enriched terms
terms <- gse$Description[1:3]
p11 <- pmcplot(terms, 2015:2024, proportion=FALSE)
p11

ggsave(filename = "figures/gse_pmcplot.pdf", plot = p11, width = 16, height = 12)



# =============================================================================
# KEGG Pathway Enrichment Analysis
# =============================================================================
print_status("\2Starting KEGG pathway enrichment analysis................line_261...........\2")
start_time <- Sys.time()
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
dedup_ids = ids[!duplicated(ids$ENSEMBL),]
df2 = df[df$ensembl_id %in% dedup_ids$ENSEMBL,]
df2$Y = dedup_ids$ENTREZID
kegg_gene_list <- df2$log2FoldChange
names(kegg_gene_list) <- df2$Y
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list <- kegg_gene_list[!duplicated(names(kegg_gene_list))]
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "hsa"
print_status("\3Running KEGG gene set enrichment analysis................line_274...........\3")
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
end_time <- Sys.time()
print_status(sprintf("KEGG analysis completed in %s", format(end_time - start_time)))

# KEGG plots
p12 <- dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
p12
ggsave(filename = "figures/KEGG_dotplot.pdf", plot = p12, width = 16, height = 12)


## Encrichment map:
### apply 'pairwise_termsim' before making plots
## calculate the similarity matrix for all gene sets, and not the top 200!
## note that this takes some time
kk3 <- pairwise_termsim(kk2, showCategory = dim(kk2)[1] )

## Extract the gene sets. These are the 20 sets that are in the dotplot (10 activated, 10 supressed).
geneSets2plot_2 <- p12$data$Description

## make a emapplot showing only the selected gene sets

p13 <- enrichplot::emapplot(kk3, showCategory=geneSets2plot_2)
p13

ggsave(filename = "figures/KEGG_emapplot.pdf", plot = p13, width = 16, height = 12)


## Category Netplot
p14 <- cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)
p14

ggsave(filename = "figures/KEGG_cnetplot.pdf", plot = p14, width = 16, height = 12)


## Ridgeplot
p15 <- ridgeplot(kk2) + labs(x = "enrichment distribution")
p15

ggsave(filename = "figures/KEGG_ridgeplot.pdf", plot = p15, width = 16, height = 12)


## GSEA Plot
# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
p16 <- gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)
p16

ggsave(filename = "figures/KEGG_gseaplot.pdf", plot = p16, width = 16, height = 12)


## KEGG Pathway Visualization with Pathview
# Produce the native KEGG plot (PNG)
# dme1 <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04130", species = "hsa")

# Produce a different plot (PDF) (not displayed here)
# dme2 <- pathview(gene.data=kegg_gene_list, pathway.id="hsa03060", species = "hsa", kegg.native = F)

# knitr::include_graphics("hsa04130.pathview.png")

# =============================================================================
# Circos Plot Visualization
# =============================================================================
print_status("\4Starting Circos plot visualization................line_357...........\4")
print_status("\5Preparing data for genome-wide expression visualization................line_358...........\5")


# Function to get human gene coordinates from Ensembl
get_human_ensembl_coordinates <- function() {
  # Connect to Ensembl
  message("Connecting to Ensembl...")
  ensembl <- useEnsembl(biomart = "genes")
  
  # Find and set human dataset
  datasets <- listDatasets(ensembl)
  human_datasets <- grep("human", datasets$description, ignore.case = TRUE)
  
  if (length(human_datasets) == 0) {
    stop("No human dataset found in Ensembl")
  }
  
  dataset <- datasets$dataset[human_datasets[1]]  # Use first human dataset
  message(paste("Using dataset:", dataset))
  
  # Set dataset and get gene coordinates
  ensembl <- useDataset(dataset = dataset, mart = ensembl)
  
  # Define attributes for gene coordinates
  attributes <- c("ensembl_gene_id", "chromosome_name", 
                  "start_position", "end_position")
  
  # Get coordinates for all genes
  message("Retrieving gene coordinates...")
  all.genes <- getBM(attributes = attributes, 
                     values = list(ensembl_gene_id = c()), 
                     mart = ensembl)
  
  # Rename column to match DESeq2 output
  colnames(all.genes)[1] <- "ensembl_id"
  
  return(all.genes)
}

# Get gene coordinates
gene_coords <- get_human_ensembl_coordinates()

# Merge coordinates with DESeq2 results
message("Merging coordinates with DESeq2 results...")
res_2 <- read.csv("results/res_2.csv")
merged_data <- merge(gene_coords, res_2, by = "ensembl_id", all.x = TRUE)

# Format chromosome names
merged_data$chromosome_name <- paste0("chr", merged_data$chromosome_name)

# Save complete dataset
output_dir <- "data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(merged_data, 
          file.path(output_dir, "deseq.csv"), 
          row.names = FALSE)

# Create and save subset of interesting genes
if (exists("genes_to_check")) {
  merged_data_subset <- merged_data[merged_data$Gene.Name %in% genes_to_check, ]
  write.csv(merged_data_subset, 
            file.path(output_dir, "deseq_subset.csv"), 
            row.names = FALSE)
  message(paste("Saved subset of", nrow(merged_data_subset), "genes"))
}

## Creating the Circos Plot
print_status("\6Creating Circos plot................line_410...........\6")
start_time <- Sys.time()


# Prepare data for Circos plot
circos_data <- merged_data %>%
  # Filter out non-standard chromosomes
  filter(grepl("^chr[0-9XY]+$", chromosome_name)) %>%
  # Sort by chromosome and position
  arrange(chromosome_name, start_position) %>%
  mutate(log2fc = ifelse(is.na(log2FoldChange), 0, log2FoldChange))

# Initialize Circos plot with chromosomes
circos.clear()
message("Initializing Circos plot...")

# Set up Circos plot parameters
circos.par(
  "track.height" = 0.1, 
  "start.degree" = 90,
  gap.degree = 2
)

# Create base chromosome ideogram
print_status("\7Initializing Circos plot with ideogram................line_438...........\7")
circos.initializeWithIdeogram(
  species = "hg19",
  chromosome.index = paste0("chr", c(1:22, "X", "Y"))
)
end_time <- Sys.time()
print_status(sprintf("Circos plot completed in %s", format(end_time - start_time)))

# Add track for log2 fold changes
circos.genomicTrack(
  circos_data[, c("chromosome_name", "start_position", "end_position", "log2fc")],
  panel.fun = function(region, value, ...) {
    circos.genomicPoints(
      region,
      value$log2fc,
      pch = 16,
      cex = 0.5,
      col = ifelse(value$log2fc > 0, "red", "blue")
    )
  },
  track.height = 0.2
)

# Add track for significant genes
message("Adding significant genes track...")
sig_genes <- circos_data %>% filter(padj < 0.05)
if (nrow(sig_genes) > 0) {
  # Prepare data for labels, ensuring correct column order
  sig_genes_for_labels <- sig_genes %>%
    dplyr::select(chromosome_name, start_position, end_position, Gene.name, log2fc)
  
  circos.genomicLabels(
    sig_genes_for_labels,
    labels.column = 4, # Gene.name column
    side = "inside",
    col = ifelse(sig_genes_for_labels$log2fc > 0, "red", "blue")
  )
}


# Save plot
pdf("figures/circos_plot.pdf", width = 10, height = 10)
dev.off()

# Print final completion message
print_status("\1Analysis pipeline completed successfully!.........line_493...........\1")
print_status(sprintf("Total execution time: %s", format(Sys.time() - as.POSIXct(Sys.getenv("START_TIME")))))
