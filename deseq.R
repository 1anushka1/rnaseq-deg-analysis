library(DESeq2)
library("ggplot2")
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
head(counts)

metadata = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/ExperimentDesignFile.RnaSeq/experiment-design")
head(metadata)

#wrangling data 
#gene id as row names
rownames(counts) = counts$Gene.ID
head(counts)

#dropping unnecessary columns
genes = counts[,c("Gene.ID", "Gene.Name")]
head(genes)
counts = counts[, -c(1,2)]
head(counts)

head(metadata)
rownames(metadata) = metadata$Run
head(metadata)

#keeping only genotype column in metadata
metadata = metadata[, "Sample.Characteristic.genotype.", drop = FALSE]
metadata

#renaming
colnames(metadata) = c("genotype")
metadata

metadata$genotype[metadata$genotype == "wild type genotype"] = "wildtype"
metadata$genotype[metadata$genotype == "Snai1 knockout"] = "knockout"
metadata

#convert to factor
metadata$genotype = factor(metadata$genotype, levels = c("wildtype", "knockout"))
metadata$genotype

colnames(counts) == rownames(metadata) #must be true


#spot check for SNAI1
#to get gene id for Snai1
gene_id = genes$Gene.ID[genes$Gene.Name == "SNAI1"]
gene_id

gene_counts = counts[gene_id,]
gene_counts

gene_data = cbind(metadata, counts = as.numeric(gene_counts))
gene_data
ggplot(gene_data, aes(x = genotype, y = counts, fill = genotype)) + geom_boxplot()

# Running DESeq
dds = DESeqDataSetFromMatrix(countData = counts, colData = metadata, design =~ genotype)

dds = dds[rowSums(counts(dds)) > 10,]
dds = DESeq(dds)

res = results(dds, contrast = c("genotype", "knockout", "wildtype"), alpha = 1e-5)
res

#making it a data frame
res_df = as.data.frame(res)
head(genes)
res_df = merge(res_df, genes, by="row.names")

#checking genes
genes_to_check = c("THY1", "SFMBT2", "PASD1", "SNAI1")
res_df[res_df$Gene.Name %in% genes_to_check,]

#MA plot
plotMA(res)

#volcano plot
library(EnhancedVolcano)

EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='pvalue')




#GO/Pathway analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

deg = res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]
genenames = deg$Gene.Name

#converting to entrez ids
genes_ensembl = deg$Gene.ID
gene_df = bitr(genes_ensembl,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

#core steps
ego = enrichGO(gene = gene_df$ENTREZID,OrgDb = org.Hs.eg.db,ont = "BP",pAdjustMethod = "BH",pvalueCutoff = 0.05)
dotplot(ego, showCategory = 10)

kegg = enrichKEGG(gene = gene_df$ENTREZID,organism = "hsa",pvalueCutoff = 0.05)
dotplot(kegg, showCategory = 10)





