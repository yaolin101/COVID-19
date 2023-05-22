rm(list = ls())
setwd(dir = "/Users/yaolin/projects_eriba/Covid_19")
options(stringsAsFactors = FALSE)
options (future.globals.maxSize = 4000 * 1024^5)

# load libraries for analysis
library(DESeq2)
library(tximeta)
library(pheatmap)
library('pvclust')
library('bitops')
library('sva')
library('limma')
library(RColorBrewer)
library(fields)
library(ggpubr)
library(tidyverse)
library(tibble)


library("AnnotationDbi")
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library(gplots)
library(dplyr)

source_colors <- c('#009E73', '#F0E442', '#E69F00', '#D55E00')

##########################SASP Mass Spectro visulization#############
library(xlsx)
SASP_meta <- read.xlsx("20210106_metadata_COVID19.xlsx", sheetIndex=1)
SASP_meta <- SASP_meta %>% 
  tidyr::separate(FGA.CODE, sep = 1, into = c("A", "B")) %>%
  dplyr::arrange(as.numeric(B))
SASP_meta$SENICAgroup <- as.factor(SASP_meta$SENICAgroup)
ggplot(SASP_meta, aes(x = SENICAgroup, y = production_age)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2))

SASP_MassSpect <- read.xlsx("DR_DEMARIA_FINAL_REPORT_JN(12-23-2021).xlsx", sheetIndex = 3, startRow = 9)
SASP_MassSpect <- SASP_MassSpect[4:48,-c(1,23,24)] # remove the useless lines and columns
SASP_MassSpect$Severity <- SASP_meta$SENICAgroup
SASP_MassSpect <- SASP_MassSpect[is.na(SASP_meta$Comments),] # remove the bad quality samples by meta data
names(SASP_MassSpect)
SASP_MassSpect <- SASP_MassSpect[,-c(3,6,14,28)] # remove the low abundance proteins

SASP_MassSpect <- lapply(SASP_MassSpect, gsub, pattern='[><]', replacement='')
SASP_MassSpect <- lapply(SASP_MassSpect, gsub, pattern='\\*', replacement='')
SASP_MassSpect <- do.call(cbind.data.frame, SASP_MassSpect)

SASP_MassSpect <- na.omit(SASP_MassSpect)
SASP_meta <- SASP_meta %>%
  filter(B %in% substring(SASP_MassSpect$RUN, 2, nchar(SASP_MassSpect$RUN))) %>%
  group_by(SENICAgroup) %>%
  summarise_at(c("production_age"), median, na.rm = TRUE)

ggplot(SASP_meta, aes(x = SENICAgroup, y = production_age)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2))
ggplot(SASP_meta, aes(x = SENICAgroup, y = gender)) +
  geom_jitter(shape=16, position=position_jitter(0.2))

annot <- as.data.frame(SASP_MassSpect$Severity); rownames(annot) <- SASP_MassSpect[,1]; colnames(annot) <- "Severity"
annoCol<-list(Severity=c(CONTROL="#009E73", NONSEVERE="#F0E442", SEVERE="#E69F00"))
source_colors <- c('#009E73', '#F0E442', '#E69F00', '#D55E00', '#56B4E9', '#000000', '#CC79A7', '#999999')

rownames(SASP_MassSpect) <- SASP_MassSpect$RUN
SASP_MassSpect <- SASP_MassSpect[,-c(1,25)]
ncol(SASP_MassSpect)
i <- 1:23
SASP_MassSpect[, i] <- apply(SASP_MassSpect[ , i], 2,           
                             function(x) as.numeric(as.character(x)))
names(SASP_MassSpect) <- c("MIF", "IGFBP1", "IGFBP3", "GDF15", "TNFRI", "Fas/CD95", "CCL2",
                           "CXCL1", "CXCL10", "CXCL2", "IL-b", "IL2", "PDGF-AA", "TNFSF10",
                           "IL-1a", "PAI1", "TNFRII", "IFNg", "IL6", "IL8", "TNF", "IL7", "IGFBP2")

# length(breaks) == length(paletteLength) + 1; use floor and ceiling to deal with even/odd length pallettelengths
paletteLength <- 200
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
dist_cor <- function(x) as.dist(1 - cor(t(x)))
max(scale(dist_cor(SASP_MassSpect))); min(scale(dist_cor(SASP_MassSpect)));
myBreaks <- c(seq(-4, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(4/paletteLength, 4, length.out=floor(paletteLength/2)))
pheatmap::pheatmap(SASP_MassSpect,scale = "column", clustering_method = "ward.D",
                   annotation_row = annot, annotation_colors = annoCol,
                   clustering_distance_rows = "correlation",cluster_cols = F,
                   color=myColor, breaks=myBreaks, show_rownames = F,
                   cutree_rows = 2, main = "SASP_MassSpect")



########################## Bulk RNA-Seq datasets ##########################
#============ load in data========
coldata <- read.csv(file = "meta_covid.csv", row.names = 1, header = T)
files <- file.path("quants_genecode_decoy",rownames(coldata), "quant.sf")
all(file.exists(files))
coldata <- data.frame(files, 
                      names = sapply(strsplit(rownames(coldata), "_"), function(x) x[2]), 
                      group = c(rep(c("control"), 10), rep(c("non_severe"), 10), rep(c("severe"), 10)))
#-------check the se object
se <- tximeta(coldata)
suppressPackageStartupMessages(library(SummarizedExperiment))
colData(se)
assayNames(se)
rowRanges(se)

#-----summatize to gene level
gse <- summarizeToGene(se)
rowRanges(gse)

gse <- addIds(gse, "SYMBOL", gene=TRUE, multiVals="list") # add different identifiers, eg RefSeq, Symbol
mcols(gse)
gse_names <- as.data.frame(mcols(gse))
#===========quality control========
dds <- DESeqDataSet(gse, ~ group)

d <- NULL
i<- 1
for (i in 1:dim(colData(gse))[1]) {
  d <- rbind(d, data.frame(x = rownames(colData(gse))[i],
                           y = log2(gse@assays@data@listData[["counts"]][, i] + 1),
                           sp = colData(gse)$group[i]))
}

yl <- expression(log[2](TPM+1))
ggplot(data = d, aes(x = x, y = y, fill = sp)) + 
  geom_boxplot(notch = TRUE, outlier.colour = "red", outlier.shape = 1) +
  labs(x = NULL, y = as.expression(yl)) +
  theme(axis.text.x = element_text(size = 8, angle = 30, hjust = 1))
dds <- dds[,-c(9,15,19,24)] # samples with extremely low sequencing depth
rld <- rlog(dds, blind=TRUE)
plotPCA(rld, intgroup="group")  + geom_label(aes(label = name))


dds <- dds[,-17] # remove the 020 sample (non-severe) because it is clustered with the control group
rld <- rlog(dds, blind=TRUE)
pcaData <- plotPCA(rld, intgroup="group", returnData = TRUE) 
percentVar <- round(100 * attr(pcaData, "percentVar")) 
ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  scale_color_manual(values = c('#009E73', '#F0E442', '#E69F00')) + 
  stat_ellipse(linetype = 2) +
  theme_bw(base_size = 15) +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#=========== senescencent markers and SASP gene expression profile ===
#normalized_counts <- counts(dds, normalized=TRUE)
#differential expression analysis
dds <- DESeq(dds)
resultsNames(dds)

#-----ENSG00000124762.14 == CDKN1A; ENSG00000147889.18 == CDKN2A
d <- plotCounts(dds, gene="ENSG00000124762.14", intgroup="group", 
                returnData=TRUE)
my_comparisons <- list(c("control", "non_severe"), c("control", "severe"), c("non_severe", "severe"))
ggplot(d, aes(x=group, y=count, color=group)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) + 
  scale_colour_manual(values=setNames(c('#009E73', '#F0E442', '#E69F00'), levels(d$group))) +
  ggtitle("CDKN1A") + 
  xlab("Severity") +
  ylab("Normalized Counts") +
  stat_compare_means(comparisons = my_comparisons)+
  theme_bw(base_size = 15) +
  theme( panel.grid.major = element_blank(), panel.border = element_blank(),
         panel.grid.minor = element_blank(), axis.text=element_text(size=15),
         axis.line = element_line(colour = "black"))

 

SASP <- c("IL1A", "IL1B", "IL6", "CXCL1", "CXCL10", "MMP1", 
          "IGFBP2", "IGFBP3", "TNF", "TIMP1", "TIMP2", "TIMP3", "TIMP4")

index <- gse_names$SYMBOL %in% SASP
mat  <- assay(rld)[ index, ]
gse_names[index,]

rownames(mat) <- gse_names[index,]$SYMBOL
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, "group"]); names(anno) <- "group"
rownames(anno) <- colnames(rld)
annoCol<-list(group=c(control="#009E73", non_severe="#F0E442", severe="#E69F00"))

pheatmap(mat, cluster_cols = F, annotation_col = anno)
pheatmap::pheatmap(mat, scale = "row", annotation_col = anno, color=myColor, breaks=myBreaks,
                   cluster_cols = F, main = "SASP gene expression", annotation_colors = annoCol)


##======================== differential expression analysis ============
res_ns <- results(dds, contrast=c("group","non_severe","control"))
res_ns_LFC <- lfcShrink(dds, contrast=c("group","non_severe","control"), type="ashr")
summary(res_ns_LFC)
write.csv(res_ns_LFC, "res_ns_LFC.csv")
resSig_ns <- res_ns_LFC %>%
  data.frame() %>%
  arrange(log2FoldChange) %>%
  dplyr::filter(padj < 0.05, abs(log2FoldChange) > 0.6)
dim(resSig_ns)
write.csv(resSig_ns, "resSig_ns_LFC.csv")

res_s <- results(dds, contrast=c("group","severe","control"))
res_s_LFC <- lfcShrink(dds, contrast=c("group","severe","control"), type="ashr")
summary(res_s_LFC)
write.csv(res_s_LFC, "res_s_LFC.csv")
resSig_s <- res_s_LFC %>%
  data.frame() %>%
  arrange(log2FoldChange) %>%
  dplyr::filter(padj < 0.05, abs(log2FoldChange) > 0.6)
dim(resSig_s)
write.csv(resSig_s, "resSig_s_LFC.csv")

res_s_ns <- results(dds, contrast=c("group","severe","non_severe"))
res_s_ns_LFC <- lfcShrink(dds, contrast=c("group","severe","non_severe"), type="ashr")
summary(res_s_ns_LFC)
write.csv(res_s_ns_LFC, "res_s_ns_LFC.csv")
resSig_s_ns <- res_s_ns_LFC %>%
  data.frame() %>%
  arrange(log2FoldChange) %>%
  dplyr::filter(padj < 0.05, abs(log2FoldChange) > 0.6)
dim(resSig_s_ns)



library("IHW")
resIHW_s <- results(dds, contrast=c("group","severe","control"), filterFun=ihw)
resIHW_ns <- results(dds, contrast=c("group","non_severe","control"), filterFun=ihw)
summary(resIHW_s)
summary(resIHW_ns)
sum(resIHW_s$padj < 0.1, na.rm=TRUE)

#------ venn diagram
library(VennDiagram)

venn.diagram(
  x = list(rownames(resSig_ns),rownames(resSig_s)),
  category.names = c("non_s" , "severe"),
  filename = 'covid_19.png'
)

lists <- list(
  DEG_s = sapply(strsplit(rownames(resSig_s),"\\."), function(x) x[1]),
  DEG_ns = sapply(strsplit(rownames(resSig_ns),"\\."), function(x) x[1])
)


merged_sig <- merge(resSig_s, resSig_ns, by="row.names", suffixes = c(".s",".ns"))

common <- get.venn.partitions(lists)[1, '..values..']
specific_ns <- get.venn.partitions(lists)[2, '..values..']
specific_s <- get.venn.partitions(lists)[3, '..values..']


##=========senescence related genes expression=========
annoCol<-list(group=c(control="#009E73", non_severe="#F0E442", severe="#E69F00"))
senmayo <- read.csv("SenMayo_human.txt", header = F)
senmayo <- unlist(senmayo$V1)


rld_mat <- as.data.frame(assay(rld))
rld_mat$Symbol <- sapply(mcols(gse, use.names=F)$SYMBOL, "[[",1)

rld_senmayo <- rld_mat[rld_mat$Symbol %in% senmayo,]
rownames(rld_senmayo) <- rld_senmayo$Symbol
rld_senmayo <- rld_senmayo[,-ncol(rld_senmayo)]
rld_senmayo <- rld_senmayo[apply(rld_senmayo, 1, function(x) !all(x==0)),]
rld_senmayo <- rld_senmayo - rowMeans(rld_senmayo)

max(rld_senmayo);min(rld_senmayo)
anno <- as.data.frame(colData(rld)[,"group", drop=FALSE])
paletteLength <- 200
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
#dist_cor <- function(x) as.dist(1 - cor(t(x)))
myBreaks <- c(seq(-4, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(4/paletteLength, 4, length.out=floor(paletteLength/2)))
pheatmap(rld_senmayo, annotation_col = anno, clustering_method = "ward.D",show_colnames=F,
         color=myColor, breaks=myBreaks,fontsize_row=5,cutree_cols=2, 
         annotation_colors = list(group=c(control="#009E73", non_severe="#F0E442", severe="#E69F00")))

##============== function analysis===============

library("org.Hs.eg.db")
library("clusterProfiler")
library(enrichplot)
library(ggnewscale)
library(GOSemSim)
library(DOSE)

# GSEA: -log10pvalue * sign(log2FoldChange)----------
library(clusterProfiler)
library(GOSemSim)
library(enrichplot)
library(msigdbr)
msigdbr_t2g <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>%
  dplyr::select(gs_name, ensembl_gene) %>% as.data.frame()
msigdbr_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, ensembl_gene) %>% as.data.frame()
msigdbr_t2g <- msigdbr(species = "Homo sapiens") %>%
  dplyr::select(gs_name, ensembl_gene) %>% as.data.frame()

# modify for each tissue and gender
matric <- -log10(res_s_LFC$pvalue) * sign(res_s_LFC$log2FoldChange)
names(matric) <- rownames(res_s_LFC)
geneList <- as.data.frame(matric) %>% 
  mutate(r = rank(matric, ties.method = "first")) %>%
  na.omit()
geneList <- geneList %>% arrange(desc(r))
geneList <- apply(geneList[,"matric", drop=FALSE], 1,c)
names(geneList) <- sapply(strsplit(names(geneList), "\\."), function(x) x[1])

em <- GSEA(geneList, TERM2GENE = msigdbr_t2g, pvalueCutoff = 0.05)
em@result$ID
em_ <- gseGO(geneList, 
             OrgDb         = org.Hs.eg.db,
             keyType       = "ENSEMBL",
             ont           = "ALL",
             pAdjustMethod = "BH")
res <- as.data.frame(em)
dotplot(em, showCategory=30)
emapplot(em_2, )
em_2 <- pairwise_termsim(em)
em_3 <- simplify(em_2)
treeplot(em_2, showCategory=30) + ggtitle("dotplot for GSEA_MG_MV_li")

write.csv(as.data.frame(em), file = paste0("GSEA_sign_MG_MV_li", ".csv" ))



#comparecluster: over representation analysis ==========
lists <- list(
  DEG_s = sapply(strsplit(rownames(subset(resSig_s,padj<0.05 && abs(log2FoldChange)>0.6)),"\\."), function(x) x[1]),
  DEG_ns = sapply(strsplit(rownames(subset(resSig_ns,padj<0.05 && abs(log2FoldChange)>0.6)),"\\."), function(x) x[1])
)
library(clusterProfiler)
library(enrichplot)
ck <- compareCluster(geneCluster = lists, 
                     OrgDb         = org.Hs.eg.db,
                     keyType="ENSEMBL",
                     fun = enrichGO,
                     ont="BP")

hsGO <- godata('org.Hs.eg.db', ont="BP")
ck_2 <- pairwise_termsim(ck)
ck_3 <- simplify(ck_2)
#emapplot(ck_3, showCategory = 20)
dotplot(ck_2, showCategory=25, font.size = 8.5) + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=40))+
  xlab("")


treeplot(ck_3)

write.csv(ck_2, "compareCluster_enrichGO.csv")






######################## modified plotPCA #############
## Modified plotPCA from DESeq2 package. Shows the Names of the Samples (the first col of SampleTable), and uses ggrepel pkg to plot them conveniently.
# @SA 10.02.2017 
library(genefilter)
library(ggplot2)
library(ggrepel)
plotPCA.san <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group, 
                  intgroup.df, name = colData(rld)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[2,3]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 
  
}
