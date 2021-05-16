setwd("/home/bioinformatica/Documentos/Martiela/TCGA_BRCA")
library("TCGAbiolinks")

query <- GDCquery(project = c("TCGA-BRCA"),
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = "Illumina HiSeq",
                  file.type  = "results", #resultados de expressao genica
                  experimental.strategy = "RNA-Seq",
                  legacy = TRUE) #HG19

# metodo api =  download em blocos (baixa de onde parou)
GDCdownload(query, method = "api", directory = "/home/bioinformatica/Documentos/Martiela/TCGA_BRCA")

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
data <- GDCprepare(query, directory = "/home/bioinformatica/Documentos/Martiela/TCGA_BRCA")

library(SummarizedExperiment)
rowRanges(data)
# cria matriz a partir de: data
dataMatrix <- assay(data,"raw_count")
# remove outliers
data_CorOutliers <- TCGAanalyze_Preprocessing(data)

# Downstream analysis using gene expression data  
# TCGA samples from IlluminaHiSeq_RNASeqV2 with type rsem.genes.results
save(data, geneInfo , file = "dataGeneExpression.rda")

# normalization of genes using dataMatrix
dataNorm <- TCGAanalyze_Normalization(tabDF = dataMatrix, geneInfo =  geneInfo)

# quantile filter of genes using normalized data
# general expression: dofferential and non differential expressed
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut =  0.25)
# normal tissue sample
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("NT"))

# tumoral tissue sample
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("TP"))

# tecido com metastase
samplesTM <- TCGAquery_SampleTypes(barcode = colnames(dataFilt),
                                   typesample = c("TM"))

# Diff.expr.analysis (DEA)
# matN: dados filtrados nas amostras NT e TP
# fdr: 0.01 - false discovery rate
# FC: 0 = sem diferenças 1 = upregulated -1 = downregulated

# Esse nao deve encontrar diferenças
# usa pra o volcano plot
# 001 = allgenes
dataDEGs_001_allgenes <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                                         mat2 = dataFilt[,samplesTP],
                                         Cond1type = "Normal",
                                         Cond2type = "Tumor",
                                         fdr.cut = 2 ,
                                         logFC.cut = 0,
                                         method = "glmLRT")

# Esse sim! E MUITA!
# Usa para enriquecimento
# 005 = 001
dataDEGs_005 <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                                mat2 = dataFilt[,samplesTP],
                                Cond1type = "Normal",
                                Cond2type = "Tumor",
                                fdr.cut = 0.01 ,
                                logFC.cut = 2,
                                method = "glmLRT")

#volcano 1
pdf("volcanoplot_001_allgenes_fc2.pdf")
with(dataDEGs_001, plot(logFC, -log10(FDR), pch=20, cex=0.5, main="Azul LogFC+TP -- Vermelho LogFC-NT", cex.main=0.5))
with(subset(dataDEGs_001, ((logFC>2) & (FDR <0.01))), points(logFC, -log10(FDR), pch=20, cex=0.3, col="blue"))
with(subset(dataDEGs_001, ((logFC<(-2)) & (FDR <0.01))), points(logFC, -log10(FDR), pch=20, cex=0.3, col="red"))
abline(v=2,col="gray60")
abline(v=-2,col="gray60")
abline(h=(-1*(log10(0.01))),col="purple")
dev.off()

#volcano
pdf("volcanoplot_005_fc2.pdf")
with(dataDEGs_all_genes, plot(logFC, -log10(FDR), pch=20, cex=0.5, main="Azul LogFC+TP -- Vermelho LogFC-NT", cex.main=0.5))
with(subset(dataDEGs_all_genes, ((logFC>2) & (FDR <0.05))), points(logFC, -log10(FDR), pch=20, cex=0.3, col="blue"))
with(subset(dataDEGs_all_genes, ((logFC<(-2)) & (FDR <0.05))), points(logFC, -log10(FDR), pch=20, cex=0.3, col="red"))
abline(v=2,col="gray60")
abline(v=-2,col="gray60")
abline(h=(-1*(log10(0.05))),col="purple")
dev.off()

# RESERVA ;D
dataDEGs_all_genes <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                                      mat2 = dataFilt[,samplesTP],
                                      Cond1type = "Normal",
                                      Cond2type = "Tumor",
                                      fdr.cut = 0.999 ,
                                      logFC.cut = 0,
                                      method = "glmLRT")

write.table(dataDEGs_001, 'dataDEGs_001_fc_1.txt', sep='\t')
write.table(dataDEGs_005, 'dataDEGs_005fc_2.txt', sep='\t')
write.table(dataDEGs_all_genes, 'dataDEGs_all_genes.txt', sep='\t')

############################################
#       ANALISE DE ENRIQUECIMENTO          #
############################################

# Usa a tabela dos diferencialmente expressos
# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(dataDEGs_005,"Tumor","Normal",
                                          dataFilt[,samplesTP],dataFilt[,samplesNT])

# Enrichment Analysis EA
# Gene Ontology (GO) and Pathway enrichment by DEGs list
Genelist <- rownames(dataDEGsFiltLevel)

system.time(ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",Genelist))

# Enrichment Analysis EA (TCGAVisualize)
# Gene Ontology (GO) and Pathway enrichment barPlot
TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = Genelist,
                        nBar = 10)

###PCA plot
pdf('pca_top_50_miRNA.pdf')
TCGAvisualize_PCA(dataFilt,dataDEGsFiltLevel,ntopgenes = 50,colnames(dataFilt[,samplesTP]),colnames(dataFilt[,samplesNT]))
###azul ficou TP e vermelho NT como no volcano
dev.off()

# Como deu erro precisa excluir s TM:
dataFilt2 <- dataFilt[,!(colnames(dataFilt)%in%samplesTM)]

###PCA plot -> apos remoção do TM
pdf('pca_top_50.pdf')
TCGAvisualize_PCA(dataFilt2,dataDEGsFiltLevel,ntopgenes = 50,colnames(dataFilt2[,samplesTP]),colnames(dataFilt2[,samplesNT]))
###azul ficou TP e vermelho NT como no volcano
dev.off()

###PCA plot -> apos remoção do TM
pdf('pca_top_100.pdf')
TCGAvisualize_PCA(dataFilt2,dataDEGsFiltLevel,ntopgenes = 100,colnames(dataFilt2[,samplesTP]),colnames(dataFilt2[,samplesNT]))
###azul ficou TP e vermelho NT como no volcano
dev.off()

write.table(dataFilt2,"samplesNT_TP.txt", sep="\t")
write.table(dataDEGs_001,"difexp_allgenes.txt", sep="\t")