#### BRCA analysis - Machine Learning posterior project #######
setwd()
library(TCGAbiolinks)
library(SummarizedExperiment)
library(EDASeq)
#To prepare gene expression data legacy = TRUE is about hg19 data

query <- GDCquery(project = c("TCGA-BRCA"),
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq",
                  file.type  = "results", 
                  experimental.strategy = "RNA-Seq",
                  legacy = TRUE) 

# api method to better download the results
GDCdownload(query, method = "api", directory = "D:/TCGA-BRCA")

#################################################################
####### Downloading data for project TCGA-BRCA ################
## GDCdownload will download 1215 files. A total of 1.842027909 GB ###
## The total size of files is big. We will download files in chunks ##
## Downloading chunk 1 of 2 (659 files, size = 999.146358 MB) as Wed_May_19_10_30_49_2021_0.tar.gz ##
###################################################################################

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
data <- GDCprepare(query, directory = "D:/TCGA-BRCA")

### To look the GDCprepare results
rowRanges(data)
#To create the data matrix with the raw counts
dataMatrix <- assay(data,"raw_count")
# remove outliers to improve the analysis
data_CorOutliers <- TCGAanalyze_Preprocessing(data)


# Downstream analysis using gene expression data  
# TCGA samples from IlluminaHiSeq_RNASeqV2 with type rsem.genes.results
save(data, geneInfo , file = "dataGeneExpression.rda")

# normalization of genes using dataMatrix
dataNorm <- TCGAanalyze_Normalization(tabDF = dataMatrix, geneInfo =  geneInfo)

# I Need about  302 seconds for this Complete Normalization Upper Quantile 
## [Processing 80k elements /s]  
## Step 1 of 4: newSeqExpressionSet ...
## Step 2 of 4: withinLaneNormalization ...
## Step 3 of 4: betweenLaneNormalization ...
## Step 4 of 4: exprs ...

# quantile filter of genes using normalized data
# general expression: differential and non differential expressed
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut =  0.25)
# normal tissue sample
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))

# tumoral tissue sample
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("TP"))

# metastatic tissue
samplesTM <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("TM"))

# Diff.expr.analysis (DEA)
# matN: filtering by samples  NT and  TP
# fdr: 0.01 - false discovery rate
# FC: 0 = sem diferenças 1 = upregulated -1 = downregulated


# To use in the volcano plot
# 001 = allgenes
dataDEGs_001_allgenes <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                                         mat2 = dataFilt[,samplesTP],
                                         Cond1type = "Normal",
                                         Cond2type = "Tumor",
                                         fdr.cut = 2 ,
                                         logFC.cut = 0,
                                         method = "glmLRT")

######################################################################
## Batch correction skipped since no factors provided               ##
## ----------------------- DEA -------------------------------      ##
##   there are Cond1 type Normal in  113 samples                    ##
## there are Cond2 type Tumor in  1095 samples                      ##
## there are  14893 features as miRNA or genes                      ##
## I Need about  600 seconds for this DEA.                          ##
##     [Processing 30k elements /s]                                 ##
######################################################################

# To find more degs
# This we need to use in the enrichment analysis
dataDEGs_005 <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                                mat2 = dataFilt[,samplesTP],
                                Cond1type = "Normal",
                                Cond2type = "Tumor",
                                fdr.cut = 0.01 ,
                                logFC.cut = 2,
                                method = "glmLRT")

#volcano plot
pdf("volcanoplot_001_allgenes_fc2.pdf")
with(dataDEGs_001_allgenes, plot(logFC, -log10(FDR), pch=20, cex=0.5, main="Blue LogFC+TP -- Red LogFC-NT", cex.main=0.5))
with(subset(dataDEGs_001_allgenes, ((logFC>2) & (FDR <0.01))), points(logFC, -log10(FDR), pch=20, cex=0.3, col="blue"))
with(subset(dataDEGs_001_allgenes, ((logFC<(-2)) & (FDR <0.01))), points(logFC, -log10(FDR), pch=20, cex=0.3, col="red"))
abline(v=2,col="gray60")
abline(v=-2,col="gray60")
abline(h=(-1*(log10(0.01))),col="purple")
dev.off()

#volcano
pdf("volcanoplot_005_fc2.pdf")
with(dataDEGs_005, plot(logFC, -log10(FDR), pch=20, cex=0.5, main="Blue LogFC+TP -- Red LogFC-NT", cex.main=0.5))
with(subset(dataDEGs_005, ((logFC>2) & (FDR <0.05))), points(logFC, -log10(FDR), pch=20, cex=0.3, col="blue"))
with(subset(dataDEGs_005, ((logFC<(-2)) & (FDR <0.05))), points(logFC, -log10(FDR), pch=20, cex=0.3, col="red"))
abline(v=2,col="gray60")
abline(v=-2,col="gray60")
abline(h=(-1*(log10(0.05))),col="purple")
dev.off()

# To explore other possible results XD
dataDEGs_all_genes <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                                      mat2 = dataFilt[,samplesTP],
                                      Cond1type = "Normal",
                                      Cond2type = "Tumor",
                                      fdr.cut = 0.999 ,
                                      logFC.cut = 0,
                                      method = "glmLRT")

# Save the expression data, filter by FoldChange
write.table(dataDEGs_001, 'dataDEGs_001_fc_1.txt', sep='\t')
write.table(dataDEGs_005, 'dataDEGs_005fc_2.txt', sep='\t')
write.table(dataDEGs_all_genes, 'dataDEGs_all_genes.txt', sep='\t')

############################################
#       ENRICHMENT ANALYSIS                #
############################################

# Here we need to provide the expression table
# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs_005,"Tumor","Normal",
                                          dataFilt[,samplesTP],dataFilt[,samplesNT])

# Enrichment Analysis EA
# Gene Ontology (GO) and Pathway enrichment by DEGs list
Genelist <- rownames(dataDEGsFiltLevel)

system.time(ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist))

#### [1] "I need about  1 minute to finish complete  Enrichment analysis GO[BP,MF,CC] and Pathways... "

# Enrichment Analysis EA (TCGAVisualize)
# Gene Ontology (GO) and Pathway enrichment barPlot
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = Genelist,
                        nBar = 20)



###blue is  tumor TP and red normal tissue 
# Maybe exclude the metastatic samples? TM:
dataFilt2 <- dataFilt[,!(colnames(dataFilt)%in%samplesTM)]

###PCA plot -> if remove metastatic samples
pdf('pca_top_50.pdf')
TCGAvisualize_PCA(dataFilt2,dataDEGsFiltLevel,ntopgenes = 50,colnames(dataFilt2[,samplesTP]),colnames(dataFilt2[,samplesNT]))
dev.off()

### top 20
pdf('pca_top_20.pdf')
TCGAvisualize_PCA(dataFilt2,dataDEGsFiltLevel,ntopgenes = 20,colnames(dataFilt2[,samplesTP]),colnames(dataFilt2[,samplesNT]))
dev.off()

### top 10
pdf('pca_top_10.pdf')
TCGAvisualize_PCA(dataFilt2,dataDEGsFiltLevel,ntopgenes = 10,colnames(dataFilt2[,samplesTP]),colnames(dataFilt2[,samplesNT]))
dev.off()

###PCA plot -> with the top 100 degs
pdf('pca_top_100.pdf')
TCGAvisualize_PCA(dataFilt2,dataDEGsFiltLevel,ntopgenes = 100,colnames(dataFilt2[,samplesTP]),colnames(dataFilt2[,samplesNT]))
dev.off()

#Save the results - again and again...
write.table(dataFilt2,"samplesNT_TP.txt", sep="\t")
write.table(dataDEGs_001_allgenes,"difexp_allgenes.txt", sep="\t")

# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs_001_allgenes,"Tumor","Normal",
                                          dataFilt[,samplesTP],dataFilt[,samplesNT])

## try to win this challenge
TCGAvisualize_Heatmap(
  dataFilt2,
  col.metadata,
  row.metadata,
  col.colors = NULL,
  row.colors = NULL,
  show_column_names = FALSE,
  show_row_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  sortCol,
  extremes = NULL,
  rownames.size = 12,
  title = NULL,
  color.levels = NULL,
  values.label = NULL,
  filename = "heatmap.pdf",
  width = 10,
  height = 10,
  type = "expression",
  scale = "none",
  heatmap.legend.color.bar = "continuous"
)


## challenge 2
TCGAVisualize_volcano(x = dataDEGs_005$logFC,
                      y = dataDEGs_005$FDR,
                      filename = "brca_volcanoexp1.png",
                      x.cut = 5,
                      y.cut = 10^-7,
                      names = rownames(dataDEGs_005),
                      color = c("black","red","darkgreen"),
                      names.size = 2,
                      xlab = " Gene expression fold change (Log2)",
                      legend = "State",
                      title = "Volcano plot (up vs down)",
                      width = 10)

### BRCA subtype information from: doi.org/10.1016/j.ccell.2018.03.014

