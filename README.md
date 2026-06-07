# RNA-seq Differential Gene Expression and Functional Enrichment Analysis

Reproduction of a published RNA-seq study investigating transcriptional changes associated with SNAI1 knockout in human cells.

## Overview
SNAI1 is a transcription factor that drives epithelial-mesenchymal transition (EMT). This analysis uses publicly available RNA-seq data (ArrayExpress: E-MTAB-5244) to identify differentially expressed genes between wild-type and SNAI1 knockout conditions, followed by pathway enrichment analysis.

## Methodology
- **Data source**: Raw counts and metadata from EBI ArrayExpress (E-MTAB-5244)
- **Differential expression**: DESeq2 (padj < 0.05, |log2FC| > 1)
- **Visualization**: MA plot, Volcano plot (EnhancedVolcano)
- **Functional enrichment**: Gene Ontology (Biological Process) and KEGG pathway analysis via ClusterProfiler

## Tools & Libraries
R, DESeq2, ClusterProfiler, EnhancedVolcano, org.Hs.eg.db, enrichplot, ggplot2
