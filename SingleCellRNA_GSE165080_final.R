# for covid project
setwd(dir = "/Users/yaolin/projects_eriba/Covid_19")
# .libPaths("/data/p303402/R/x86_64-pc-linux-gnu-library/4.2")

install.packages(c("BiocManager","Seurat","SCINA","ggplot2","devtools","shiny"))
BiocManager::install(c("SingleCellExperiment","scater","dplyr",
                       "scmap","celldex","SingleR"))
devtools::install_github("romanhaa/cerebroApp")
devtools::install_github("mw201608/msigdb")


library(dplyr)
library(Seurat)
require(openxlsx)
require(ggplot2)
library(ggpubr)
require(cowplot)
library(data.table)
library(RColorBrewer)
library(SingleR)
library(scater)
library(pheatmap)
library(tidyverse)
library(rstatix)

#-------------GSE165080------
# 53 patient totally including healthy control
GSE165080.counts <- Read10X("/Users/yaolin/projects_eriba/Covid_19/paper/GSE165080")
GSE165080 <- CreateSeuratObject(counts = GSE165080.counts,
                                min.cells = 3,
                                min.genes = 300,
                                min.features = 300,
                                names.delim = "-"
)

rm(GSE165080.counts) # to save spaces
#VlnPlot(object = GSE165080, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), cols = 3,  pt.size = 0.1)
#FeatureScatter(object = GSE165080, "nCount_RNA", "nFeature_RNA")
#FeatureScatter(object = GSE165080, "percent.mito", "nFeature_RNA")


phynotype <- read.csv("/Users/yaolin/projects_eriba/Covid_19/paper/GSE165080/phynotype.tsv", 
                      sep = "", header = F, row.names = 1)

cell_batch <- read.csv("/Users/yaolin/projects_eriba/Covid_19/paper/GSE165080/GSE165080_cell_batch.tsv",
                       sep = "\t")

phynotype <- as.data.frame(t(phynotype))


GSE165080$subject <- plyr::mapvalues(colnames(GSE165080), from = cell_batch$Cell, to=cell_batch$batch)
GSE165080$subject_age <- plyr::mapvalues(GSE165080$subject, from = phynotype$Sample, to = phynotype$Sample_age)
GSE165080$subject_gender <- plyr::mapvalues(GSE165080$subject, from = phynotype$Sample, to = phynotype$Sample_gender)
GSE165080$subject_severity <- plyr::mapvalues(GSE165080$subject, from = phynotype$Sample, to = phynotype$Sample_severity)
GSE165080$subject_duration <- plyr::mapvalues(GSE165080$subject, from = phynotype$Sample, to = phynotype$Sample_duration)
GSE165080$subject_comorbidity <- plyr::mapvalues(GSE165080$subject, from = phynotype$Sample, to = phynotype$Sample_comorbidity)


GSE165080_subset <- subset(GSE165080, 
                           subset =  (subject_severity == "Health_control" | subject_severity ==  "Moderate" | subject_severity ==  "Severe"))
GSE165080_subset2 <- subset(GSE165080_subset, 
                            subset = (subject_comorbidity == "No" | subject_comorbidity == "Health_control"))

rm(GSE165080_subset)
rm(GSE165080)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = GSE165080_subset2), value = TRUE)
percent.mito <- colSums(GSE165080_subset2[mito.genes,]) / colSums(GSE165080_subset2)
GSE165080_subset2 <- AddMetaData(object = GSE165080_subset2, metadata = percent.mito, col.name = "percent.mito")

GSE165080_subset3 <- subset(GSE165080_subset2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito < 0.25)

#SCTransform as an alternative to the NormalizeData, FindVariableFeatures, ScaleData workflow.
#Transformed data will be available in the SCT assay, which is set as the default after running sctransform
#During normalization, we can also remove confounding sources of variation, for example, mitochondrial mapping percentage
#glmGamPoi package substantially improves the speed of the learning procedure. 
GSE165080_subset2 <- SCTransform(GSE165080_subset2, method = "glmGamPoi", vars.to.regress = "percent.mito", conserve.memory = T, verbose = T)

GSE165080_subset2_test <- Seurat::NormalizeData(GSE165080_subset2) # Normalize the data
GSE165080_subset2_test <- Seurat::FindVariableFeatures(GSE165080_subset2_test) # Determine the variable features of the dataset
GSE165080_subset2_test <- Seurat::ScaleData(GSE165080_subset2_test) # Scale the data based on the variable features


#seurat <- FilterCells(object = seurat, subset.names = c("nGene", "percent.mito","n.exp.hkgenes"), low.thresholds = c(350, -Inf,55), high.thresholds = c(5000, 0.1, Inf))
#Perform dimensionality reduction by PCA and UMAP embedding
GSE165080_subset2_test <- RunPCA(GSE165080_subset2_test, verbose = FALSE)
GSE165080_subset2_test <- RunUMAP(GSE165080_subset2_test, dims = 1:30, verbose = FALSE)

GSE165080_subset2_test <- FindNeighbors(GSE165080_subset2_test, dims = 1:30, verbose = FALSE)
GSE165080_subset2_test <- FindClusters(GSE165080_subset2_test, verbose = FALSE)
DimPlot(GSE165080_subset2_test, label = TRUE) + NoLegend()

#================================== annotate the cell types by SingleR =====
# SingleR() expects reference datasets to be normalized and log-transformed.
# For the test data, the assay data need not be log-transformed or even (scale) normalized. 
set.seed(9742)
ref <- celldex::MonacoImmuneData()

#query_sce <- Seurat::as.SingleCellExperiment(GSE165080_subset2_test)
query_sce <- Seurat::as.SingleCellExperiment(GSE165080_subset3)
query_sce <- scater::logNormCounts(query_sce)
predictions <- SingleR::SingleR(test=query_sce, ref=ref, labels=ref$label.fine)

# Add singleR labels to query_sce
sum(is.na(predictions$pruned.labels))
predictions$pruned.labels[which(is.na(predictions$pruned.labels))] <- "ambiguous"
colData(query_sce)$singleR <- predictions$pruned.labels

# Create a bar plot of number of cells per assigned cell ID
par(mar=c(13, 4, 2, 0))
table(predictions$pruned.labels)
s
# Make a UMAP
query_sce <- scater::runUMAP(query_sce)

save(query_sce, file="query_sce_from_SingleCellRNA_GSE165080.RData")

#============================== visulization ======================

load("/Users/yaolin/projects_eriba/Covid_19/query_sce_from_SingleCellRNA_GSE165080.RData")

query.seurat <- as.Seurat(query_sce, counts = "counts", data = "logcounts")
query.seurat <- ScaleData(query.seurat, features = rownames(query.seurat))
query.seurat$subject_severity <- gsub("Health_control", "Health", query.seurat$subject_severity)

Idents(query.seurat) <- "singleR"
DimPlot(query.seurat, label = TRUE, repel = T) + NoLegend()

query.seurat.mono <- subset(query.seurat, idents=c("Intermediate monocytes", "Classical monocytes", "Non classical monocytes"))
#query.seurat.ncmono <- subset(query.seurat.mono, idents=c("Intermediate monocytes", "Non classical monocytes"))

#Idents(query.seurat) <- "subject_severity"
barplot(table(query.seurat$singleR), las=2, angle = 45); test <- table(query.seurat$singleR)

plot <- DimPlot(query.seurat)
LabelClusters(plot = plot, id = "ident")

FeaturePlot(query.seurat, features = "CDKN1A")
FeaturePlot(query.seurat, features = c("CDKN1A", "CD14", "CD16", "LYZ"))
VlnPlot(query.seurat.mono, features = "CDKN1A", split.by = "subject_severity", pt.size = 0, 
        cols = c('#009E73', '#F0E442', '#E69F00')) +
  theme(axis.text.x=element_text(angle = 0,hjust=0.5, size=15)) + xlab("")

################################# cell cycle position #################
library(tricycle)

query_sce_mono <- query_sce[,query_sce$singleR %in% c("Intermediate monocytes", "Classical monocytes", "Non classical monocytes")]
query_sce_mono$subject_severity <- gsub("Health_control", "Health", query_sce_mono$subject_severity)

query_sce_mono <- project_cycle_space(query_sce_mono, species = "human", gname.type = "SYMBOL")
query_sce_mono
scater::plotReducedDim(query_sce_mono, dimred = "tricycleEmbedding") +
  labs(x = "Projected PC1", y = "Projected PC2") +
  ggtitle(sprintf("Projected cell cycle space (n=%d)",
                  ncol(query_sce_mono))) +
  theme_bw(base_size = 14)

query_sce_mono <- estimate_cycle_position(query_sce_mono)


#cell cycle density
names(colData(query_sce_mono))
colors <- c('#009E73', '#F0E442', '#E69F00')
library(Cairo)
cairo_pdf(file = "~/Desktop/test.pdf", width = 8, height = 5)
plot_ccposition_den(query_sce_mono$tricyclePosition,line.size=1,line.alpha=0.8,
                    factor(query_sce_mono$subject_severity), 'severity',
                    palette.v = c('#009E73', '#F0E442', '#E69F00'),
                    bw = 10, fig.title = "Density of Cell Cycle position") +
  theme_bw(base_size = 14)+
  theme( panel.grid.major = element_blank(),legend.justification = "top",
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

cairo_pdf(file = "~/Desktop/test2.pdf", width = 8, height = 5)
circle_scale_legend(text.size = 5, alpha = 0.9, addStageLabel = T)
dev.off()



#============ AUCell ===============
library(AUCell)
library(ggsignif)


senmayo <- read.csv("SenMayo_human.txt", header = F)
senmayo <- unlist(senmayo$V1)
geneSets <- list(SenMayo=senmayo)

exprMatrix <- GetAssayData(query.seurat.mono, slot = "data") # slot = "counts"
cells_AUC <- AUCell_run(exprMatrix, geneSets)
#cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=TRUE, splitByBlocks=TRUE)
#cells_AUC <- AUCell_calcAUC(senmayo, cells_rankings)

query.seurat.mono <- AddMetaData(query.seurat.mono, t(as.data.frame(cells_AUC@assays@data@listData[["AUC"]])), col.name = "AUC" )
violin_data <- query.seurat.mono[[c("subject_severity", "singleR","AUC")]]

ggplot(violin_data, aes(factor(singleR),AUC, fill="subject_severity")) +
  geom_violin(trim=TRUE,aes(fill=factor(subject_severity)),show.legend = T, position = position_dodge(0.9)) + 
  geom_boxplot(aes(fill=factor(subject_severity)),width = 0.15, outlier.shape = NA, position = position_dodge(0.9) ) + 
  scale_fill_manual(values = c('#009E73', '#F0E442', '#E69F00'))+
  theme_bw() + xlab("")+ ylab("SenMayo AUC score") +
  theme(axis.text.x=element_text(size=15))+ 
  theme_bw(base_size = 15) +
  theme(legend.title=element_blank())+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# add p value
my_comparisons <- list(c("Healthy_control", "Moderate"), c("Healthy_control", "Severe"))

anno_df = compare_means(AUC ~ subject_severity, group.by = "singleR", data = violin_data) %>%
  mutate(y_pos = 0.2)

ggplot(violin_data, aes(factor(subject_severity),AUC, fill="subject_severity")) +
  geom_violin(trim=TRUE,aes(fill=factor(subject_severity)),show.legend = T, position = position_dodge(0.9) ) + 
  geom_boxplot(aes(fill=factor(subject_severity)),width = 0.15, outlier.shape = NA, position = position_dodge(0.9) ) + 
  facet_wrap(~singleR) + 
  ggsignif::geom_signif(
    data=anno_df, 
    aes(xmin=group1, xmax=group2, annotations=p.adj, y_position=y_pos), 
    manual=TRUE) +
  theme_bw() + xlab("") +
  theme_bw(base_size = 15) +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# another method
stat.test <- violin_data %>%
  group_by(singleR) %>%
  t_test(AUC ~ subject_severity, ref.group = "Health_control") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "singleR", dodge = 0.8)


ggviolin(violin_data, x = "singleR", y = "AUC", color = "subject_severity", palette = c('#009E73', '#F0E442', '#E69F00'))  + 
  stat_pvalue_manual(stat.test,  label = "p", tip.length = 0, hide.ns = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))



# p21 potive and negative separately, to see if the SenMayo are from the p21 positive cells.
my_comparisons_p21 <-  list(c("CDKN1A_NEG", "CDKN1A_POS"))

poscells <- WhichCells(query.seurat.mono, expression = CDKN1A > 0)
query.seurat.mono$CDKN1A_status <- ifelse(colnames(query.seurat.mono) %in% poscells, "CDKN1A_POS", "CDKN1A_NEG")

violin_p21_data <- query.seurat.mono[[c("subject_severity", "singleR", "AUC", "CDKN1A_status")]]
violin_p21_data <- violin_p21_data[violin_p21_data$subject_severity %in% "Severe",]

anno_df_p21 = compare_means(AUC ~ CDKN1A_status, group.by = "singleR", data = violin_p21_data) %>%
  mutate(y_pos = 0.2)

ggplot(violin_p21_data, aes(x=CDKN1A_status, y=AUC)) + 
  geom_violin(trim=TRUE,aes(fill=factor(CDKN1A_status)),show.legend = T, position = position_dodge(0.9) ) + 
  geom_boxplot(aes(fill=factor(CDKN1A_status)),width = 0.15, outlier.shape = NA, position = position_dodge(0.9) ) + 
  scale_fill_manual(values = c('#009E73', '#F0E442'))+
  theme_bw() + ylab("SenMayo AUC Score") +
  facet_wrap(~singleR) + 
  ggsignif::geom_signif(
    data=anno_df_p21, 
    aes(xmin=group1, xmax=group2, annotations=p.adj, y_position=y_pos), 
    manual=TRUE) +
  theme_bw(base_size = 15) +
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


#save.image(file='singleCellRNA_GSE165080_v2.RData')
#load('singleCellRNA_GSE165080_v2.RData')

