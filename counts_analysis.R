# Analyze the neurospora crassa expression data for a collaboration with emily

# import the libraries we need for the analysis

library(dplyr)
library(DESeq2)
library(tximport)
library(readr)
library(GenomicFeatures)
library(dplyr)
library(stringr)
library(ggrepel)
library(tidyverse)
library(pheatmap)
library(DESeq2)
library(apeglm)
library(RColorBrewer)
 




# Read in the annotation data

coldata <- data.table::fread("Annotation/Sample_Info.csv", header = TRUE)
rownames(coldata) <- coldata$UUID # check colnames in read in


# DONE
# read in the counts data and make an unstranded RNA counts
samples <- list.files("~/storage/neurospora_crassa/geneCounts")
i <- samples[1]
tmp <- dat[,1]
tmp <- tmp[,1]
for(i in samples){
  dat <- data.table::fread(stringr::str_glue("~/geneCounts/{i}"), header = FALSE)
  colnames(dat) <- c("geneID","unstranded","forwardStrand","reverseStrand")
  # cut f1,f2
  
  dat <- dat %>%
    dplyr::select(unstranded)
  
  colnames(dat) <- str_sub(i,1,16)
  
  tmp <- cbind(tmp,dat)
  
  
}
#write.csv(x = df, file = "/mnt/storage/mskaro1/neurospora_crassa/results/mRNA_counts.csv")


######################################
######################################
######################################
######################################

# Analyze differential expression of NC genes
df <- df[5:length(df$geneID),]
colnames(df) <- c("geneID","4420_29CON_A_S21","4420_29CON_B_S22","4420_29CON_C_S23","4420_29CON_D_S24","4420_29MYC_A_S9R","4420_29MYC_B_S10","4420_29MYC_C_S11","4420_29MYC_D_S12","4420_72CON_A_S13","4420_72CON_B_S14","4420_72CON_C_S15","4420_72CON_D_S16","4420_72MYC_A_S1","4420_72MYC_B_S2","4420_72MYC_C_S3","4420_72MYC_D_S4","4420_76CON_A_S17","4420_76CON_B_S18","4420_76CON_C_S19","4420_76CON_D_S20","4420_76MYC_A_S5","4420_76MYC_B_S6","4420_76MYC_C_S7","4420_76MYC_D_S8")

rownames(coldata) <- coldata$UUID
coldata$Strain <- as.factor(coldata$Strain)
coldata$Cell_Type <- as.factor(coldata$Cell_Type)

# assert the gene ID as the rownames

df <- df %>%
  tibble::column_to_rownames("geneID")

coldata <- coldata %>%
  dplyr::select(c("Strain", "Cell_Type"))

rownames(coldata) %in% colnames(df)
rownames(coldata) == colnames(df)

dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = coldata,
                              design = ~ Strain + Cell_Type)
featureData <- data.frame(gene=rownames(df))
mcols(dds) <- DataFrame(mcols(dds), featureData)
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]
#dds$condition <- factor(dds$condition, levels = c("Tumor","Normal"))
#dds$condition <- relevel(dds$condition, ref = "Normal")

dds <- DESeq(dds)

res <- results(dds)

resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")

for(i in 2:length(resultsNames(dds))){
  #print(i)
  j <- resultsNames(dds)[i]
  resLFC <- lfcShrink(dds, coef=resultsNames(dds)[i], type="apeglm")
  
  resOrdered <- as.data.frame(resLFC[order(resLFC$pvalue),])
  #write.csv(x = resOrdered, file = str_glue("{j}_differential_expression_analysis.csv"))
  
  # DGE
  
  de <- resOrdered[complete.cases(resOrdered), ]
  
  # Convert directly in the aes()
  # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
  
  # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
  
  # add a column of NAs
  de$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  de$diffexpressed[de$log2FoldChange > 5 & de$padj < 1e-25] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  de$diffexpressed[de$log2FoldChange < -5 & de$padj < 1e-25] <- "DOWN"
  
  
  de$delabel <- NA
  de$delabel[de$padj<1e-25 & abs(de$log2FoldChange) > 5] <- rownames(de)[de$padj<1e-25 & abs(de$log2FoldChange) > 5]
  
  p <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
    geom_point() + 
    theme_minimal() +
    geom_vline(xintercept=c(-5, 5), col="red") +
    geom_hline(yintercept=-log10(10e-27), col="red") + 
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=c("blue", "black", "red"))
  
  print(j)
  
  ggsave(filename = str_glue("Significantly_differenially_expressed_genes_using_power_FC_for_fdr_{j}.pdf"),plot = p,device = "pdf",width = 8,
         height = 6,units = "in")
  
}


save(file = "mRNA_DGE.Rdata", dds)
write.csv(model.matrix(design(dds), colData(dds)),file = "model_mat.csv")


coldata <- data.table::fread("~/storage/neurospora_crassa/Annotation/Sample_Info.csv", header = TRUE)
rownames(coldata) <- coldata$UUID # check colnames in read in


# DONE
# read in the counts data and make an unstranded RNA counts
samples <- list.files("~/storage/neurospora_crassa/geneCounts")
i <- samples[1]
tmp <- dat[,1]
tmp <- tmp[,1]
for(i in samples){
  dat <- data.table::fread(stringr::str_glue("~/storage/neurospora_crassa/geneCounts/{i}"), header = FALSE)
  colnames(dat) <- c("geneID","unstranded","forwardStrand","reverseStrand")
  # cut f1,f2
  
  dat <- dat %>%
    dplyr::select(unstranded)
  
  colnames(dat) <- str_sub(i,1,16)
  
  tmp <- cbind(tmp,dat)
  
  
}
#write.csv(x = df, file = "/mnt/storage/mskaro1/neurospora_crassa/results/mRNA_counts.csv")


######################################
######################################
######################################
######################################
df <- data.table::fread("/mnt/storage/mskaro1/neurospora_crassa/results/mRNA_counts.csv") %>%
  dplyr::select(-V1) %>%
  tibble::column_to_rownames("geneID")
# Analyze differential expression of NC genes
df <- df[5:length(df[,1]),]
colnames(df) <- c("4420_29CON_A_S21","4420_29CON_B_S22","4420_29CON_C_S23","4420_29CON_D_S24","4420_29MYC_A_S9R","4420_29MYC_B_S10","4420_29MYC_C_S11","4420_29MYC_D_S12","4420_72CON_A_S13","4420_72CON_B_S14","4420_72CON_C_S15","4420_72CON_D_S16","4420_72MYC_A_S1","4420_72MYC_B_S2","4420_72MYC_C_S3","4420_72MYC_D_S4","4420_76CON_A_S17","4420_76CON_B_S18","4420_76CON_C_S19","4420_76CON_D_S20","4420_76MYC_A_S5","4420_76MYC_B_S6","4420_76MYC_C_S7","4420_76MYC_D_S8")

rownames(coldata) <- coldata$UUID
coldata$Strain <- as.factor(coldata$Strain)
coldata$Cell_Type <- as.factor(coldata$Cell_Type)

# assert the gene ID as the rownames


# strain and cell type
cell_type_strain <- function(coldata, df){
  

coldata <- coldata %>%
  dplyr::select(c("Strain", "Cell_Type"))

rownames(coldata) %in% colnames(df)
rownames(coldata) == colnames(df)

dds <- DESeqDataSetFromMatrix(countData = df,
                              colData = coldata,
                              design = ~ Strain + Cell_Type)
featureData <- data.frame(gene=rownames(df))
mcols(dds) <- DataFrame(mcols(dds), featureData)
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]
#dds$condition <- factor(dds$condition, levels = c("Tumor","Normal"))
#dds$condition <- relevel(dds$condition, ref = "Normal")

dds <- DESeq(dds)

res <- results(dds)

resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")

for(i in 2:length(resultsNames(dds))){
  #print(i)
  j <- resultsNames(dds)[i]
  resLFC <- lfcShrink(dds, coef=resultsNames(dds)[i], type="apeglm")
  
  resOrdered <- as.data.frame(resLFC[order(resLFC$pvalue),])
  #write.csv(x = resOrdered, file = str_glue("{j}_differential_expression_analysis.csv"))
  
  # DGE
  
  de <- resOrdered[complete.cases(resOrdered), ]
  
  # Convert directly in the aes()
  # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
  
  # The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
  # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
  
  # add a column of NAs
  de$diffexpressed <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  de$diffexpressed[de$log2FoldChange > 5 & de$padj < 1e-25] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  de$diffexpressed[de$log2FoldChange < -5 & de$padj < 1e-25] <- "DOWN"
  
  
  de$delabel <- NA
  de$delabel[de$padj<1e-25 & abs(de$log2FoldChange) > 5] <- rownames(de)[de$padj<1e-25 & abs(de$log2FoldChange) > 5]
  
  p <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
    geom_point() + 
    theme_minimal() +
    geom_vline(xintercept=c(-5, 5), col="red") +
    geom_hline(yintercept=-log10(10e-27), col="red") + 
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=c("blue", "black", "red"))
  
  print(j)
  
  ggsave(filename = str_glue("Significantly_differenially_expressed_genes_using_power_FC_for_fdr_{j}.pdf"),plot = p,device = "pdf",width = 8,
         height = 6,units = "in")
  
}


save(file = "mRNA_DGE.Rdata", dds)
write.csv(model.matrix(design(dds), colData(dds)),file = "model_mat.csv")


# make visualizations of each of the comparison tables

vsd <- vst(dds, blind=FALSE)

head(assay(vsd), 3)

# distance matrix
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Strain, vsd$Cell_Type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("Strain", "Cell_Type"))

pcaData <- plotPCA(vsd, intgroup=c("Strain", "Cell_Type"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Strain, shape=Cell_Type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(values=c("#E69F00", "#68838B", "#CD6839"))


#Strain 8872: lightblue4 (#68838B)
#Strain 8876: sienna3 (#CD6839)
#Strain 2229: navajowhite (#FFDEAD)# this is a terribe color, sorry have to change


library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Strain", "Cell_Type")])
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df, show_colnames = FALSE)

mcols(dds)$devSat <- mcols(dds)$deviance - 2*rowSums(dnbinom(counts(dds), mu=counts(dds), size=1/dispersions(dds), log=TRUE))
mc = as.data.frame(mcols(dds))
write.csv(x = mc, file = "All_pairwise_comprisons_statistical_analysis_sheet.csv")
# OK, the issue is that -2 times the log likelihood needs the -2LL(saturated model) 
# subtracted out before it can be used as a quality of fit statistic (if at all). 
# Which not hard to do with the existing output. If I might make a request/suggestion, 
# it would be great to add an additional column to the output of DESeq (the mcols-accessed slot) 
# that has -2LL(fitted model) -2LL(saturated model). Perhaps the new quantity could be called devianceVsSaturatedModel.

plot(log(mc$baseMean),log(mc$deviance))
plot(log(mc$baseVar),log(mc$deviance))



colnames(All_pairwise_comprisons_statistical_analysis_sheet)

}


strain <- function(coldata, df){
  coldata <- coldata %>%
    dplyr::filter(coldata$Cell_Type == "Conidiophore") %>%
    dplyr::select("Strain")
  
  df <- df %>%
    dplyr::select(rownames(coldata))
  
  rownames(coldata) == colnames(df)
  
  dds <- DESeqDataSetFromMatrix(countData = df,
                                colData = coldata,
                                design = ~ Strain)
  featureData <- data.frame(gene=rownames(df))
  mcols(dds) <- DataFrame(mcols(dds), featureData)
  keep <- rowSums(counts(dds)) >=10
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- results(dds)
  resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
  resOrdered <- as.data.frame(resLFC[order(resLFC$pvalue),])
  
  save(file = "strain_DGE_conidiophore.Rdata", dds)
  write.csv(model.matrix(design(dds), colData(dds)),file = "Conidiophore_model_mat_Strain.csv")
  write.csv(x = resOrdered, file = "Conidiophore_Strain_8872_vs_2229_res.csv")
  
  
  resLFC <- lfcShrink(dds, coef=resultsNames(dds)[3], type="apeglm")
  resOrdered <- as.data.frame(resLFC[order(resLFC$pvalue),])
  
  write.csv(x = resOrdered, file = "Conidiophore_Strain_8876_vs_2229_res.csv")
  
  mcols(dds)$devSat <- mcols(dds)$deviance - 2*rowSums(dnbinom(counts(dds), mu=counts(dds), size=1/dispersions(dds), log=TRUE))
  mc = as.data.frame(mcols(dds))
  write.csv(x = mc, file = "All_pairwise_comprisons_Conidiophore_Strain.csv")
  
  vsd <- vst(dds, blind=FALSE)
  
  pcaData <- plotPCA(vsd, intgroup=c("Strain"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=Strain)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    scale_color_manual(values=c("#E69F00", "#68838B", "#CD6839"))
  
  library("pheatmap")
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:20]
  df <- as.data.frame(colData(dds)["Strain"])
  pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
           cluster_cols=FALSE, annotation_col=df, show_colnames = FALSE)
  
}



# make visualizations of each of the comparison tables

vsd <- vst(dds, blind=FALSE)

head(assay(vsd), 3)

# distance matrix
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Strain, vsd$Cell_Type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("Strain", "Cell_Type"))

pcaData <- plotPCA(vsd, intgroup=c("Strain", "Cell_Type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Strain, shape=Cell_Type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(values=c("#E69F00", "#68838B", "#CD6839"))


#Strain 8872: lightblue4 (#68838B)
#Strain 8876: sienna3 (#CD6839)
#Strain 2229: navajowhite (#FFDEAD)# this is a terribe color, sorry have to change


library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Strain", "Cell_Type")])
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df, show_colnames = FALSE)

mcols(dds)$devSat <- mcols(dds)$deviance - 2*rowSums(dnbinom(counts(dds), mu=counts(dds), size=1/dispersions(dds), log=TRUE))
mc = as.data.frame(mcols(dds))
write.csv(x = mc, file = "All_pairwise_comprisons_statistical_analysis_sheet.csv")
# OK, the issue is that -2 times the log likelihood needs the -2LL(saturated model) 
# subtracted out before it can be used as a quality of fit statistic (if at all). 
# Which not hard to do with the existing output. If I might make a request/suggestion, 
# it would be great to add an additional column to the output of DESeq (the mcols-accessed slot) 
# that has -2LL(fitted model) -2LL(saturated model). Perhaps the new quantity could be called devianceVsSaturatedModel.

plot(log(mc$baseMean),log(mc$deviance))
plot(log(mc$baseVar),log(mc$deviance))




##########GSEA
library(clusterProfiler)
library(enrichplot)
library(AnnotationHub)
library(ggridges)
library(enrichplot)
library(DOSE)

hub <- AnnotationHub() #creating an AnnotationHub object
query(hub, "Neurospora") #auery whether there is an annotation of Neurospora in the Database
NC <- hub[["AH95468"]] #Alwyas choose the data which have eg.sqlite extension
                       #Otherwise the following steps won't work


###Creating a gene_list for gseGO and enrichGO functions
# we want the log2 fold change of all genes from res_df dataframe
original_gene_list <- res_df$log2FoldChange

# name the vector
names(original_gene_list) <- res_df$EntrezID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

gene_list <- sort(gene_list, decreasing = TRUE) # sort the list in decreasing order (required for clusterProfiler)

###Creating the gseGO object
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "ENTREZID", 
             #nPerm = 10000, 
             #minGSSize = 3, 
             #maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = FALSE, 
             OrgDb = NC, 
             pAdjustMethod = "BH") #Optimizing the parameters as per my need

#dotplot
dotplot(gse, showCategory=10, title = 'Enriched Biological Processes', split = ".sign", font.size = 10) + facet_grid(.~.sign)

#emap plot
pairwise_termsim (gse) %>%
  emapplot(showCategory = 20, cex_label_category = 0.7)

#cnet plot
cnetplot (gse, showCategory = 3)

#ridge plot
ridgeplot(gse, showCategory = 20) + 
  labs(x = "Enrichment Distribution") +
  theme (axis.text.y = element_text (size = 10, hjust = 1, color = 'steelblue')) 

#gseaplot
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)


###creating genes for creating an enrichGO object
sig_genes_df <- subset(res_df, padj < 0.05) #choosing only the significant genes from res_df
                                            #the gene_list contains info of all the genes
                                            #the genes contains info of only the significant genes
genes <- sig_genes_df$log2FoldChange

# Name the vector
names(genes) <- sig_genes_df$EntrezID

# omit NA values
genes <- na.omit(genes)

# filter on min log2fold change (log2FoldChange > 1.58)
genes <- names(genes)[abs(genes) > 1.58] #choosing  only the names of those genes whose lfc is > 1.58

###creating the enrichGO object
ego <- enrichGO(gene = genes,
                universe = names (gene_list),
                OrgDb = NC,
                keyType = 'ENTREZID',
                readable = TRUE,
                ont = 'BP')

#cnetplot
cnetplot(ego, showCategory = 3, foldChange=gene_list, circular = TRUE, colorEdge = TRUE) 



#####KEGG pathway analysis
library(pathview)

###Creating a gseKEGG object
kegg_ncr <- gseKEGG(geneList     = gene_list,
                    organism     = "ncr", #organism in ncr, for different organisms different, check kegg database
                    nPerm        = 10000,
                    minGSSize    = 3,
                    maxGSSize    = 800,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    keyType       = "ncbi-geneid") #the keyType is ncbi-geneid 

#dotplot
dotplot(kegg_ncr, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

#emapplot
pairwise_termsim (kegg_ncr) %>%
  emapplot(showCategory = 20)

#cnetplot
cnetplot(kegg_ncr)

###PathView 

#the value of pathway.id will create a pathway.png file with that name
#go to the current working directory to find the image

ncr_TCA <- pathview(gene.data = gene_list, pathway.id = 'ncr00020',     
                      species = "ncr", limit = list(gene = 2.5, cpd = 2)) #the value of gene is the limit of abs(lfc) I want to keep

#downloading only the pathway
knitr::include_graphics('ncr00020.pathview.png')








