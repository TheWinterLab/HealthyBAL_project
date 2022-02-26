source("BAL_functions.R")
source("web_browse.R")

load("20191219_bal_all_v3.Robj")
bal.all.control <- renameClusters_2(bal.all.control,
                                    cluster.names  = c("Macrophages", "Macrophages", "T/NK Cells",
                                                       "Macrophages", "Macrophages", "Macrophages",
                                                       "Macrophages", "Dendritic Cells", "Macrophages",
                                                       "Dividing Cells", "Club Cells", "B Cells", 
                                                       "Ciliated Cells"),
                                    annotation.name = "Seurat_Assignment")

# Top 3 Genes Heatmap -----------------------------------------------------

bal.all.control.labeled.markers <- FindAllMarkers(object = bal.all.control, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(bal.all.control.labeled.markers %>% group_by(cluster), file="20190917_all_am_labeled_markers.txt", row.names=FALSE, sep="\t")
bal.all.control.labeled.markers <- read.delim(file="20190917_all_am_labeled_markers.txt")
cluster.averages <- AverageExpression(object = bal.all.control, return.seurat = TRUE)
top3 <- bal.all.control.labeled.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

CairoPNG("20200721_Top3_Genes_heatmap.png", bg="transparent",
         width = 768, height = 768)
(DoHeatmap(cluster.averages, features = top3$gene, 
           draw.lines = FALSE, size = 6, group.colors = c("#00a9ff","#cd9600","#00eea1","#e4ff00","#bf5e58","#c77cff","#ff61cc")) + viridis::scale_fill_viridis()
)
dev.off()


# FeaturePlot -------------------------------------------------------------

CairoPNG("202001118_FeaturePlot_CD68_MRC1_EPCAM.png", bg="transparent",
         width = 1500, height = 750)

#(FeaturePlot(bal.all.control, features = c("FABP4","LTB","MARCKS","MS4A1","SLPI","SCGB3A1","CAPS","PCLAF"),
            # ncol = 2, cols = c("grey", "red"), coord.fixed = TRUE, pt.size = .1)
#)

(FeaturePlot(bal.all.control, features = c("FABP4", "CCL5", "LTB", "MS4A1", "CTSB", "SCGB3A1", "CAPS", "TOP2A"),
 ncol = 4, cols = c("grey", "red"), coord.fixed = FALSE, pt.size = .1)
)

dev.off()

#CD68, MRC1, CDH1, EPCAM

CairoPNG("20201118_FeaturePlot_EPCAM.png", bg="transparent",
         width = 375, height = 375)

(FeaturePlot(bal.all.control, features = c("EPCAM"),
             ncol = 1, cols = c("grey", "red"), pt.size = .1)
)
dev.off()



# Umap with cell numbers --------------------------------------------------

CairoPNG(file = "20200825_Unlabeled_UMAP_with_Cell_Numbers.png",bg = "transparent",
         width = 1500, height = 1500)

(
  umapCellAnno(bal.all.control,
               label.size = 15, point.size = 1, cell.legend.size = 12)
)

dev.off()




# Cluster Correlation -----------------------------------------------------

cluster.averages <- AverageExpression(object = bal.all.control, show.progress = TRUE)

combined.corr<-cor(cluster.averages$SCT,method="pearson")
pheatmap(combined.corr)

CairoPNG("20200721_Cluster_Correlation.png", bg="transparent",
         width = 576, height = 576)
(
  pheatmap(combined.corr)
)

dev.off()



# Unlabled UMAP with Cell Numbers -----------------------------------------

load("20191219_bal_all_v3.Robj")

CairoPNG(file = "20200721_Unlabeled_UMAP_with_Cell_Numbers.png",bg = "transparent",
         width = 960, height = 960)

(
  umapCellAnno(bal.all.control, counts.as.title = TRUE, 
               label.size = 17, point.size = .1,
               axis.text.x.bottom.size  = 40,
               axis.text.y.left.size = 40,
               axis.title.y.size = 37,
               axis.title.x.size = 37,
               legend.text.size = 20)
)

dev.off()

theme(text = element_text(size = 40),
      axis.text.x = element_text(size = 35),
      axis.text.y = element_text(size = 35)) + 
  guides(colour = guide_legend(override.aes = list(size = 7))) 



# Quality Control Figures -------------------------------------------------


bal01 <- CreateSeuratObject(counts = Read10X("20191209_BAL01_counts_GRCh38_CRV3/filtered_feature_bc_matrix"), 
                            min.cells = 3, min.features = 200, project = "bal_ctrl_01")
bal01 <- PercentageFeatureSet(bal01, pattern = "^MT-", col.name = "percent.mt")
mito.low.cutoff <- quantile(bal01$percent.mt, probs = 0.050)
mito.high.cutoff <- quantile(bal01$percent.mt, probs = 0.975)
feature.cutoff <- quantile(bal01$nFeature_RNA, probs = 0.975)

plot1 <- FeatureScatter(bal01, feature1 = "nCount_RNA", feature2 = "percent.mt",
                        cols = "#f8766d") +
  geom_hline(yintercept=mito.low.cutoff, linetype="dashed", color = "red") +
  geom_hline(yintercept=mito.high.cutoff, linetype="dashed", color = "red") +
  theme(text = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40)) + 
  ggtitle("")
plot2 <- FeatureScatter(bal01, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        cols = "#f8766d") +
  geom_hline(yintercept=feature.cutoff, linetype="dashed", color = "red")  +
  theme(text = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40)) + 
  ggtitle("")
plots.combined <- CombinePlots(plots = list(plot1, plot2), legend = "none") + draw_figure_label(label = "BAL01")

CairoPDF("20201012_bal01_qc.png", bg="transparent",
         width = 10, height = 960)
(
  plots.combined
)
dev.off()

bal01 <- subset(bal01, subset = nFeature_RNA < feature.cutoff & percent.mt > mito.low.cutoff & percent.mt < mito.high.cutoff)



bal02.data <- Read10X("20191209_BAL02_counts_GRCh38_CRV3/filtered_feature_bc_matrix")
bal02 <- CreateSeuratObject(counts = bal02.data, min.cells = 3, min.features = 200, project = "bal_ctrl_02")
bal02 <- PercentageFeatureSet(bal02, pattern = "^MT-", col.name = "percent.mt")


mito.low.cutoff <- quantile(bal02$percent.mt, probs = 0.050)
mito.high.cutoff <- quantile(bal02$percent.mt, probs = 0.975)
feature.cutoff <- quantile(bal02$nFeature_RNA, probs = 0.975)

plot1 <- FeatureScatter(bal02, feature1 = "nCount_RNA", feature2 = "percent.mt",
                        cols = "#7cae00") +
  geom_hline(yintercept=mito.low.cutoff, linetype="dashed", color = "#7cae00") +
  geom_hline(yintercept=mito.high.cutoff, linetype="dashed", color = "#7cae00") +
  theme(text = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40)) + 
  ggtitle("")
plot2 <- FeatureScatter(bal02, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        cols = "#7cae00") +
  geom_hline(yintercept=feature.cutoff, linetype="dashed", color = "#7cae00")  +
  theme(text = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40))+
  ggtitle("")
plots.combined <- CombinePlots(plots = list(plot1, plot2), legend = "none") + draw_figure_label(label = "BAL02")

CairoPDF("202001012_bal02_qc.png", bg="transparent",
         width = 10, height = 960)
(
  plots.combined
)
dev.off()

bal02 <- subset(bal02, subset = nFeature_RNA < feature.cutoff & percent.mt > mito.low.cutoff & percent.mt < mito.high.cutoff)




bal03.data <- Read10X("20191209_BAL03_counts_GRCh38_CRV3/filtered_feature_bc_matrix")
bal03 <- CreateSeuratObject(counts = bal03.data, min.cells = 3, min.features = 200, project = "bal_ctrl_03")
bal03 <- PercentageFeatureSet(bal03, pattern = "^MT-", col.name = "percent.mt")


mito.low.cutoff <- quantile(bal03$percent.mt, probs = 0.050)
mito.high.cutoff <- quantile(bal03$percent.mt, probs = 0.975)
feature.cutoff <- quantile(bal03$nFeature_RNA, probs = 0.975)

plot1 <- FeatureScatter(bal03, feature1 = "nCount_RNA", feature2 = "percent.mt",
                        cols = "#00bfc4") +
  geom_hline(yintercept=mito.low.cutoff, linetype="dashed", color = "#00bfc4") +
  geom_hline(yintercept=mito.high.cutoff, linetype="dashed", color = "#00bfc4") +
  theme(text = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40))+
  ggtitle("")
plot2 <- FeatureScatter(bal03, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        cols = "#00bfc4") +
  geom_hline(yintercept=feature.cutoff, linetype="dashed", color = "#00bfc4")  +
  theme(text = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40))+
  ggtitle("")
plots.combined <- CombinePlots(plots = list(plot1, plot2), legend = "none") + draw_figure_label(label = "BAL03")

CairoPDF("20201012_BAL03_qc.png", bg="transparent",
         width = 10, height = 960)
(
  plots.combined
)
dev.off()

bal03 <- subset(bal03, subset = nFeature_RNA < feature.cutoff & percent.mt > mito.low.cutoff & percent.mt < mito.high.cutoff)


bal04.data <- Read10X("20191209_BAL04_counts_GRCh38_CRV3/filtered_feature_bc_matrix")
bal04 <- CreateSeuratObject(counts = bal04.data, min.cells = 3, min.features = 200, project = "bal_ctrl_04")
bal04 <- PercentageFeatureSet(bal04, pattern = "^MT-", col.name = "percent.mt")


mito.low.cutoff <- quantile(bal04$percent.mt, probs = 0.050)
mito.high.cutoff <- quantile(bal04$percent.mt, probs = 0.975)
feature.cutoff <- quantile(bal04$nFeature_RNA, probs = 0.975)

plot1 <- FeatureScatter(bal04, feature1 = "nCount_RNA", feature2 = "percent.mt",
                        cols = "#c77cff") +
  geom_hline(yintercept=mito.low.cutoff, linetype="dashed", color = "#c77cff") +
  geom_hline(yintercept=mito.high.cutoff, linetype="dashed", color = "#c77cff") +
  theme(text = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40))+
  ggtitle("")
plot2 <- FeatureScatter(bal04, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                        cols = "#c77cff") +
  geom_hline(yintercept=feature.cutoff, linetype="dashed", color = "#c77cff")  +
  theme(text = element_text(size = 40),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40))+
  ggtitle("")
plots.combined <- CombinePlots(plots = list(plot1, plot2), legend = "none") + draw_figure_label(label = "bal04")

CairoPDF("20201012_BAL04_qc.png", bg="transparent",
         width = 10, height = 960)
(
  plots.combined
)
dev.off()

bal04 <- subset(bal04, subset = nFeature_RNA < feature.cutoff & percent.mt > mito.low.cutoff & percent.mt < mito.high.cutoff)




# Unlabeled UMAP ----------------------------------------------------------
load("20191219_bal_all_v3.Robj")

p1 <- DimPlot(bal.all.control, label = TRUE, reduction = "umap", label.size = 12)  +
  theme(legend.position = "none") #+
  #ggtitle("Unlabeled Unmerged Clusters")

CairoPNG("20200721_Unlabeled_UMAP.png", bg="transparent",
         width = 960, height = 960)
(
  p1
)
dev.off()



# Clustering by Participant -----------------------------------------------
CairoPNG(file = "20201014_Unlabeled_UMAP_Grouped_By_Cell_Identity.png",bg = "transparent",
         width = 1500, height = 1500)

(
  DimPlot(bal.all.control, group.by = "orig.ident") + 
    scale_color_discrete(name = "Individuals", 
                         labels = c("BAL01","BAL02", "BAL03","BAL04")) + 
    theme(text = element_text(size = 40),
          axis.text.x = element_text(size = 35),
          axis.text.y = element_text(size = 35)) + 
    guides(colour = guide_legend(override.aes = list(size = 7))) 
) 
dev.off()




# Top 3 Gene Heatmap Unlabeled --------------------------------------------

bal.all.control.labeled.markers <- FindAllMarkers(object = bal.all.control, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(bal.all.control.labeled.markers %>% group_by(cluster), file="20190917_all_am_labeled_markers.txt", row.names=FALSE, sep="\t")
bal.all.control.labeled.markers <- read.delim(file="20190917_all_am_labeled_markers.txt")
cluster.averages <- AverageExpression(object = bal.all.control, return.seurat = TRUE)
top3 <- bal.all.control.labeled.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

CairoPNG("20200721_Top3_Genes_Unlabeled_heatmap.png", bg="transparent",
         width = 768, height = 768)
(DoHeatmap(cluster.averages, features = top3$gene, 
           draw.lines = FALSE, size = 6) + viridis::scale_fill_viridis()
)
dev.off()


# Unlabeled Cluster Correlation ---------------------------------------------------

cluster.averages <- AverageExpression(object = bal.all.control, show.progress = TRUE)

combined.corr<-cor(cluster.averages$SCT,method="pearson")
pheatmap(combined.corr)

CairoPNG("20200721_Unlabeled_Cluster_Correlation.png", bg="transparent",
         width = 576, height = 576)
(
  pheatmap(combined.corr)
)

dev.off()



# Dot Plot ----------------------------------------------------------------

cell.table <- as_tibble(FetchData(bal.all.control, vars = c("orig.ident","Seurat_Assignment")))
counts <- cell.table %>% group_by(orig.ident, Seurat_Assignment) %>% tally()

counts <- counts %>% mutate(Percentage = case_when(
  orig.ident == "bal_ctrl_01" ~ n / (cell.table %>% group_by(orig.ident) %>% tally(wt = n))$n[[1]] * 100,
  orig.ident == "bal_ctrl_02" ~ n / (cell.table %>% group_by(orig.ident) %>% tally(wt = n))$n[[2]] * 100,
  orig.ident == "bal_ctrl_03" ~ n / (cell.table %>% group_by(orig.ident) %>% tally(wt = n))$n[[3]] * 100,
  orig.ident == "bal_ctrl_04" ~ n / (cell.table %>% group_by(orig.ident) %>% tally(wt = n))$n[[4]] * 100
    ))

p1 <- ggplot(count, aes(x=Seurat_Assignment, y=Percentage)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize = .5) + 
  scale_y_continuous(trans='log2', breaks = c(.1,.5,1,2,5,10,20,50,80,100),limits = c(.1,100.0)) +
  labs(x = "Cell Type") + 
  scale_fill_discrete(name = "Individual", labels = c("BAL01", "BAL02", "BAL03","BAL04")) + 
  theme(text = element_text(size = 40))  + 
  geom_boxplot(alpha = .1) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
 
CairoPNG("20200728_CellType_Percentage_by_Individual.png", bg="transparent",
         width = 1500, height = 1500)
(
 p1
)

dev.off()


# Dot Plot with Colors ----------------------------------------------------


counts %>% group_by(Seurat_Assignment) %>% 
  mutate(min=min(Percentage,na.rm=T),max=max(Percentage,na.rm=T),avg=mean(Percentage,na.rm=T)) -> df2

set.seed(5)

p2 <- ggplot(df2, aes(x=Seurat_Assignment, y=Percentage, fill = factor(orig.ident))) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 1, position = position_jitter(0.05), alpha = .5) + 
  geom_errorbar(aes(ymin = min, ymax = max), width = 0.5,size=1)+
  geom_errorbar(aes(ymin = min, ymax = avg), width = 0.5,size=1)+
  scale_fill_discrete(name = "Individual", labels = c("BAL01", "BAL02", "BAL03","BAL04")) + 
  scale_y_continuous(trans='log2', breaks = c(.05,.1,.5,1,2,5,10,20,50,80,100),limits = c(.05,100),
                     labels = scales::comma)+
  theme(text = element_text(size = 38)) +
  labs(x = "Cell Type") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(panel.background = element_rect(fill = "#FFFFFF")) +
  theme(panel.border = element_rect(fill = NA))


CairoPNG("20200805_CellType_Percentage_by_Individual_With_Color.png", bg="transparent",
         width = 1500, height = 1500)
(
  p2
)

dev.off()



# Volcano Plot ------------------------------------------------------------

#Use Code in figures.R
# load("20200729_bal_mp.Robj")
# 
# cluster.names = c("Macrophages","MIP-1 Macrophages")
# names(cluster.names) <- levels(x = bal.mp)
# bal.mp <- Seurat::RenameIdents(bal.mp, cluster.names)
# bal.mp@meta.data$Seurat_Assignment <- Idents(bal.mp)
# 
# 
# 

# bal.mp <- FindVariableFeatures(bal.mp, nfeatures = 50)
# p1 <- VariableFeaturePlot(bal.mp)
# 
# p2 <- LabelPoints(p1, points= head(VariableFeatures(bal.mp), 10), size = 8, repel = T) + 
#   theme(text = element_text(size = 20),
#         axis.ticks = element_line(size = 2),
#         axis.title.x = element_text(size = 30),
#         axis.title.y = element_text(size = 30),
#         axis.text.x = element_text(size = 20),
#         axis.text.y = element_text(size = 20),
#         panel.border = element_rect(fill = "NA"),
#         legend.position = "none")
# 
# 
# CairoPNG("20201014_Macrophage_Variable_Gene_Plot.png", bg="transparent",
#          width = 1500, height = 1500)
# (
#   p2
# )
# 
# dev.off()


# Variable Gene Matrix ----------------------------------------------------

cor.mat <- variableGeneMatrix(bal.mp, variable.genes = 50, return.table = TRUE)
cor.mat.plot <- pheatmap(cor.mat, angle_col = 90, fontsize = 13)

CairoPNG("20200729_Variable_Gene_Correlation_Matrix.png", bg="transparent",
         width = 1500, height = 1500)
(
  cor.mat.plot
)

dev.off()

# Macrophage UMAP ---------------------------------------------------------

p1 <- umapCellAnno(bal.mp)

CairoPNG(file = "20200728_Labeled_Macrophages_UMAP_with_Cell_Numbers.png",bg = "transparent",
         width = 660, height = 720)


(
  umapCellAnno(bal.mp, 
               label.size = 7, point.size = .1, cell.legend.size = 12)
)

dev.off()

# MP Module Score ---------------------------------------------------------

bal.mp <- AddModuleScore(object = bal.mp,features = list(c("CCL3", "CCL4", "CXCL10", "CCL4L2" ,"CCL20")),
                         name = "Module_Scores")

CairoPNG("20200831_Macrophage_MIP_Subcluster_Module_Score.png", bg="transparent",
         width = 375, height = 375)

(
  FeaturePlot(bal.mp, features = "Module_Scores1", pt.size = .1)
  + ggtitle("MIP-1 Module Score") + 
    theme(text = element_text(size = 15)) + 
    theme(axis.text = element_text(size = 20))
)

dev.off()


# Feature Plot for Macrophages --------------------------------------------
split.bal.obj <- SplitObject(bal.mp, split.by = "orig.ident")

features <- c("CD68", "MRC1", "CCL4", "CXCL10", "CCL20", "CCL4L2") 



CairoPNG("20201118_FeaturePlot_FAPB4.png", bg="transparent",
         width = 375, height = 375)

(FeaturePlot(bal.mp, features = c("FABP4"),
             ncol = 1, cols = c("grey", "red"), pt.size = .1)
)
dev.off()

CairoPNG("20200826_FeaturePlot_CCL3.png", bg="transparent",
         width = 375, height = 375)

(FeaturePlot(bal.mp, features = c("CCL3"),
             ncol = 1, cols = c("grey", "red"), pt.size = .1)
)
dev.off()

CairoPNG("20200826_FeaturePlot_CCL4.png", bg="transparent",
         width = 375, height = 375)

(FeaturePlot(bal.mp, features = c("CCL4"),
             ncol = 1, cols = c("grey", "red"), pt.size = .1)
)
dev.off()


# Differential Genes MIP and Macrophages/Variable Gene Heatmap ------------

mip.all.markers <- FindMarkers(bal.mp, ident.1 = "MIP-1 Macrophages", ident.2 = "Macrophages")
mip.markers <- cbind(gene = row.names(mip.all.markers), mip.all.markers)

top5.mip.markers <- mip.markers %>% top_n(n = 10, wt = avg_logFC)
bot5.mip.markers <- mip.markers %>% top_n(n = 10, wt = -avg_logFC)
am.cells <- colnames(x = subset(x = bal.mp,
                                idents = c("MIP-1 Macrophages", "Macrophages"),
                                downsample = 500))


CairoPNG("20201014_Macrophage_MIP_Var_Gene_Heatmap_tophalf.png", bg="transparent",
         width = 1500, height =1500)

#Heatmap of top differential genes
#(DoHeatmap(bal.mp,
#          features = c(as.character(top5.mip.markers$gene), as.character(bot5.mip.markers$gene)),
#          cells = am.cells) +scale_fill_viridis())

#Heatmap based upon top variable genes
(
  DoHeatmap(bal.mp, 
            cells = am.cells, 
            features = rownames(top5.mip.markers), 
            size = 12,
            angle = 0,
            hjust = .5) +
    scale_fill_viridis() +
    theme(text = element_text(size = 20))  +
    guides(fill = guide_legend(title="Relative Expression", reverse = TRUE, guides = FALSE),
           color = guide_colorbar())
) 



dev.off()

CairoPNG("20200914_Macrophage_MIP_Var_Gene_Heatmap_bothalf.png", bg="transparent",
         width = 1500, height =1500)

(
  DoHeatmap(bal.mp, 
            cells = am.cells, 
            features = rownames(bot.mip.markers), 
            angle = 0,
            hjust = .5,
            size = 10) + 
    scale_fill_viridis() +
    theme(text = element_text(size = 20))
) 

dev.off()

# Over Representation Bar Plot --------------------------------------------


all.genes <- read_tsv("20191218_mip1_gsea.rnk", col_names = c("gene_name","fold_change"))
mip.markers <- mip.marker$gene
entrez.all.genes <- mapIds(org.Hs.eg.db, all.genes$gene_name, 'ENTREZID', 'SYMBOL')
entrez.mip.markers <- mapIds(org.Hs.eg.db, mip.markers, 'ENTREZID', 'SYMBOL')


ego2 <- enrichGO(gene         = entrez.mip.markers,
                 OrgDb         = org.Hs.eg.db,
                 universe = entrez.all.genes,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)



CairoPNG("20200729_Barplot_Over_Representation_Test.png", bg="transparent",
         width =600, height = 600)

(barplot(ego2, showCategory = 20))

dev.off()



# BAL05 FeaturePlot -------------------------------------------------------

load("bal05.Robj")

CairoPNG("20200828_BAL05_FeaturePlot_CCL4.png", bg="transparent",
         width = 1500, height = 1500)

#(FeaturePlot(bal.all.control, features = c("FABP4","LTB","MARCKS","MS4A1","SLPI","SCGB3A1","CAPS","PCLAF"),
# ncol = 2, cols = c("grey", "red"), coord.fixed = TRUE, pt.size = .1)
#)

(FeaturePlot(bal05, features = c("SPP1","FN1","CD68"),
             ncol = 3, cols = c("grey", "red"), coord.fixed = TRUE, pt.size = .1)
)

dev.off()

# BAL05 Dot Plot with Colors ----------------------------------------------

cell.table <- as_tibble(FetchData(bal05, vars = c("Seurat_Assignment")))
cell.table$Seurat_Assignment <- as.character(cell.table$Seurat_Assignment)
cell.table$Seurat_Assignment <- str_replace(cell.table$Seurat_Assignment, "MIP-1 Macrophages","Macrophages")
cell.table.total <- cell.table %>% group_by(Seurat_Assignment) %>% tally() %>% mutate(Percentage = (n/sum(count$n)) * 100)



set.seed(5)

p2<- ggplot(cell.table.total, aes(x=Seurat_Assignment, y = Percentage, fill = Seurat_Assignment)) + 
  geom_bar(stat = "identity") +
  theme(panel.background = element_rect(fill = "#ffffff")) +
  theme(text = element_text(size = 30)) +
  theme(panel.border = element_rect(fill = NA)) +
  labs(x = "Cell Type") +
  scale_fill_discrete(name = "Cell Type") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

CairoPNG("202000806_CellType_Percentage_BAL05_With_Color.png", bg="transparent",
         width = 1500, height = 1500)
(
  p2
)

dev.off()





# Module Scores BAL05 -----------------------------------------------------

bal05 <- Seurat::AddModuleScore(object = bal05,
                                features = c("CCL3", "CCL4", "CXCL10", "CCL4L2" ,"CCL20"),
                                name = "Module_Scores")

CairoPNG("20200818_BAL05_Module_Score.png", bg="transparent",
         width = 1500, height = 1500)

(
  FeaturePlot(bal05, features = "Module_Scores1", pt.size = 2) + ggtitle("") + 
    theme(text = element_text(size = 20)) + 
    theme(axis.text = element_text(size = 20)) + ggtitle("MIP-1 Module")
)

dev.off()


# Ambient RNA Scatter Plot ------------------------------------------------


sample.names <- c("BAL01","BAL02","BAL03","BAL04")
sample.scatter <- map(sample.names, sampleScatter)
sample.scatter2 <- map2(sample.scatter, sample.names, ~ plot.grid(.x,.y))

CairoPNG("202000807_Ambient_RNA_Scatter_Plot.png", bg="transparent",
         width = 1500, height = 1500)
(
  cowplot::plot_grid(sample.scatter2[[1]],sample.scatter2[[2]], sample.scatter2[[3]], sample.scatter2[[4]],
                     nrow = 2, ncol =2)
) 

dev.off()

celltypeGraph()


# Ambient RNA Included/Excluded Heatmap -----------------------------------


# excluded_included_counts <- map(c("BAL01","BAL02","BAL03","BAL04"), excluded_included_heatmap)
# cor_mat <- excluded_included_counts %>% 
#   purrr::reduce(full_join, by = "gene_name") %>%
#   as_tibble() %>%
#   filter(across(.col = -gene_name, ~.x > 0)) %>%
#   dplyr::select(-gene_name) %>%
#   mutate(across(.cols = everything(), ~replace_na(.x, 0))) %>%
#   log() %>%
#   rstatix::cor_mat() %>%
#   dplyr::select(rowname, all_of(str_subset(colnames(.), pattern = "excluded"))) %>%
#   dplyr::filter(rowname %in% str_subset(rowname, pattern = "included")) %>%
#   column_to_rownames(var = "rowname") 
# 
# pheatmap(cor_mat, border_color = NA, cluster_rows = FALSE, cluster_cols = FALSE)

original.ident <- unique(bal.all.control$orig.ident)
matrix.sample.list <- extractMatrixSums(samples = c("BAL01","BAL02","BAL03","BAL04"))
hmap <- fourByFour(original.ident = original.ident, sample.matrix = matrix.sample.list)
rownames(hmap) <- c("BAL01","BAL02","BAL03","BAL04")
colnames(hmap) <- c("BAL01","BAL02","BAL03","BAL04")

(
CairoPNG("202000902_Ambient_RNA_Excluded_Included_By_Individual_Heatmap.png", bg="transparent",
         width = 1500, height = 1500)

  
) 

dev.off()


# HL Donor Macrophages Module Scores --------------------------------------

Macrophages.All.subset <- AddModuleScore(object = Macrophages.All.subset,
                                         features = list(c("CCL3", "CCL4", "CXCL10", "CCL4L2" ,"CCL20")),
                                         name = "Module_Scores")

CairoPNG("20200831_HL_Macrophage_Donor_Module_Score.png", bg="transparent",
         width = 375, height = 375)

(
  FeaturePlot(Macrophages.All.subset, features = "Module_Scores1", pt.size = .1)
  + ggtitle("MIP-1 Module Score") 
)

dev.off()



# HL Donor Macrophages Feature Plot ---------------------------------------

CairoPNG("20200826_HL_Macrophages_FeaturePlot_CCL4.png", bg="transparent",
         width = 375, height = 375)

(FeaturePlot(Macrophages.All.subset, features = c("CCL4"),
             ncol = 1, cols = c("grey", "red"), pt.size = .1)
)

dev.off()


# HL Donor Macrophages UMAP Grouped By Participant ------------
#Add Donor Numbers to Object

load("20200811_Macrophages_Donors_Only.Robj")
Macrophages.All.subset <- RunUMAP(Macrophages.All.subset, dims = 1:30, seed.use = 24)

#Adding Donor IDS to object instead of renaming Idents

subject_id <- FetchData(Macrophages.All.subset, vars = "subject.id") %>% as_tibble(rownames = "Cells")

subject_mapping <- tibble(subject.id = unique(Macrophages.All.subset$subject.id),
                          donor.id = str_c("Donor",1:8))

subject_id <- subject_id %>% full_join(subject_mapping) %>% dplyr::select(-subject.id) %>% 
  deframe(.) %>%
  factor(levels = str_c("Donor",1:8))

Macrophages.All.subset$donor_id <- subject_id

load("HL.expressors.Robj")
load("HL.non.expressors.Robj")


exp.vec <- rep("Module High",1160)
names(exp.vec) <- HL.expressors %>% pull(name) %>% unique()

nexp.vec <- rep("Module Low", 12456)
names(nexp.vec) <-  HL.non.expressors %>% pull(name) %>% unique()
Macrophages.All.subset$module_group <- c(exp.vec, nexp.vec)





p1 <- DimPlot(Macrophages.All.subset, group.by = "donor_id")

CairoPNG("20200928_HL_Macrophage_Donor_UMAP_Grouped_by_Participant.png", bg="transparent",
         width = 375, height = 375)

(
  p1
)

dev.off()


# Hierarchical Clustering ---------------------------------------------------

load("20200729_bal_mp.Robj")
bal.mp <- renameClusters(bal.mp, cluster.names = c("Macrophages","MIP-1 Macrophages"),
                         annotation.name = "Seurat_Assignment")

mac <- hier(cluster = "Macrophages", seed = 2, down.sample = 2000)
mip <- hier(cluster = "MIP-1 Macrophages", seed = 2, down.sample = 400)

p1 <- pheatmap(t(mac[[1]]), cluster_cols = mac[[2]], show_colnames = FALSE, color = viridis(10))
p2 <- pheatmap(t(mip[[1]]), cluster_cols = mip[[2]], show_colnames = FALSE, color = viridis(10))

data <- t(rbind(mac[[1]][mac[[2]]$order,],mip[[1]][mip[[2]]$order,]))
annotation <- data.frame(Group = c(rep("Macrophages",1000), rep("MIP-1 Macrophages",200)))
rownames(annotation) <- colnames(data)
p3 <- pheatmap(t(rbind(mac[[1]][mac[[2]]$order,],mip[[1]][mip[[2]]$order,])), cluster_cols = FALSE, annotation_col = annotation,
               color = viridis(10), show_colnames = FALSE)

CairoPNG("20200814_Macrophages_MIP_Dendro_Heatmap.png", bg="transparent",
         width = 1500, height = 1500)

(
  p1
) 

dev.off()


# Violin Plot .3 Module Score Expressors ----------------------------------
#This section was run on Quest due to RAM limitations

# library(Seurat)
# library(tidyverse)
# 
# load("20200811_Macrophages_Donors_Only.Robj")
# Macrophages.All.subset <- RunUMAP(Macrophages.All.subset, dims = 1:30, seed.use = 24)
# Macrophages.All.subset <- AddModuleScore(object = Macrophages.All.subset,
#                                          features = list(c("CCL3", "CCL4", "CXCL10", "CCL4L2" ,"CCL20")),
#                                          name = "Module_Scores")
# 
# mod.scores.HL <- as_tibble(FetchData(Macrophages.All.subset, vars = "Module_Scores1"), rownames = "Cells")
# HL.expressors <- (mod.scores.HL %>% filter(Module_Scores1 >= .3))$Cells %>%
#   subset(Macrophages.All.subset, cells = .) %>%
#   GetAssayData() %>%
#   as_tibble(rownames = "Genes") %>%
#   pivot_longer(cols = -Genes) %>%
#   filter(Genes == "CCL3" | Genes == "FABP4") %>%
#   mutate(Group = "Expressors")
# 
# save(HL.expressors, file = "HL.expressors.Robj")
# 
# 
# HL.non.expressors <- (mod.scores.HL %>% filter(Module_Scores1 < .3))$Cells %>%
#   subset(Macrophages.All.subset, cells = .) %>%
#   GetAssayData() %>%
#   as_tibble(rownames = "Genes") %>%
#   pivot_longer(cols = -Genes) %>%
#   filter(Genes == "CCL3" | Genes == "FABP4") %>%
#   mutate(Group = "Non-Expressors")

# save(HL.non.expressors, file = "HL.non.expressors.Robj")

load("HL.expressors.Robj")
load("HL.non.expressors.Robj")

HL.expressors.and.non <- bind_rows(HL.expressors,HL.non.expressors) %>%
  mutate(Group = factor(Group, levels = c("Non-Expressors","Expressors")))


# CairoPNG("20200904_Expressors_Non_Expressors_ViolPlot_FABP4.png", bg="transparent",
#          width = 1500, height = 1500)
# 
# (
# ggplot(data = filter(HL.expressors.and.non, Genes == "FABP4"), mapping = aes(x = Genes, y = value, fill = Group)) +
#   geom_violin() + 
#   geom_jitter(size = 1) +
#   theme(text = element_text(size = 40),
#         panel.border = element_rect(fill = "NA"),
#         panel.background = element_rect(fill = "NA")) + 
#   labs(y = "Expression") +
#   theme(legend.position = "none") + 
#   facet_wrap(~Group, scales = "free") 
# ) 
# 
# dev.off()
# 
# CairoPNG("20200904_Expressors_Non_Expressors_ViolPlot_CCL3.png", bg="transparent",
#          width = 1500, height = 1500)
# 
# (
# ggplot(data = filter(HL.expressors.and.non, Genes == "CCL3"), 
#              mapping = aes(x = Genes, y = value, fill = Group)) +
#   geom_violin() + 
#   geom_jitter(size = 1) +
#   theme(text = element_text(size = 40),
#         panel.border = element_rect(fill = "NA"),
#         panel.background = element_rect(fill = "NA")) + 
#   facet_wrap(~Group, scales = "free") +
#     theme(legend.position = "none") + 
#   labs(y = "Expression")
# ) 
# 
# dev.off


# HL Donor UMAP by Participant AND High/Low Module -----------------------------------------

load("HL.expressors.Robj")
load("HL.non.expressors.Robj")

HL.expressors.and.non <- bind_rows(HL.expressors,HL.non.expressors) %>%
  mutate(Group = factor(Group, levels = c("Non-Expressors","Expressors")))


load("20200811_Macrophages_Donors_Only.Robj")
Macrophages.All.subset <- RunUMAP(Macrophages.All.subset, dims = 1:30, seed.use = 24)

exp.vec <- rep("Module High",1160)
names(exp.vec) <- HL.expressors %>% pull(name) %>% unique()

nexp.vec <- rep("Module Low", 12456)
names(nexp.vec) <-  HL.non.expressors %>% pull(name) %>% unique()
Macrophages.All.subset$module_group <- c(exp.vec, nexp.vec)

donor_split <- SplitObject(Macrophages.All.subset, "donor_id")

exp_count <- map(donor_split, ~ FetchData(.x, "module_group") %>% 
  as_tibble(rownames = "ID") %>% 
  group_by(module_group) %>% 
  tally() %>%
  pull(n) %>%
  paste0(c("Module High: ", "Module Low: "), .))

exp_count <- unname(exp_count)

for (i in seq_along(exp_count)){ 
  exp_count[[i]][[2]] <- paste0(exp_count[[i]][[2]], "\n")
  }

exp_count <- map_chr(exp_count, ~ paste0(.x, collapse = " \n "))
exp_count <- paste0(paste0(paste0("Donor", 1:8), "\n "), exp_count)
part_umap <- DimPlot(Macrophages.All.subset, group.by = "donor_id", pt.size = 1) + 
                     theme(text = element_text(size = 20),
                           axis.text = element_text(size  = 20))

g <- ggplot_build(part_umap )
values <- unique(g$data[[1]]["colour"])[,1]

part_umap + scale_color_manual(values = values, labels = exp_count)

CairoPDF("20201118_HL_Donor_UMAP_By_Donor_and_Module.pdf", bg="transparent",
         width = 13, height = 10)

(
  part_umap + scale_color_manual(values = values, labels = exp_count)
)

dev.off()

CairoPNG("20201118_HL_Donor_UMAP_By_Donor_and_Module.png",
         width = 1200, height = 1000)

(
  part_umap + scale_color_manual(values = values, labels = exp_count)
)

dev.off()



enhanced_vln_plot <- function(seurat.obj, features){ 
 p1 <- VlnPlot(seurat.obj, group.by = "module_group", features = features, cols = c("#cfcbd6","#6f45f6")) + 
    theme(text = element_text(size = 40)) +
    theme(axis.text.x = element_text(size = 40)) +
    theme(axis.text.y = element_text(size = 40)) + 
    theme(legend.position = "none") + 
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank())
  current.date <- str_replace_all(Sys.Date(), "-","")
  file.name <- glue("{current.date}_{features}_Expressors_Non_Expressors_VlnPlot.png")

  CairoPNG(file.name, bg="transparent",
           width = 1500, height = 1500)

  ggsave(file.name, plot = p1, width = 20, height = 20)
}

# 
# 
# walk(c("CD68","FABP4", "CCL3", "CCL4", "CCL4L2", "CCL20", "CXCL10"), ~ 
#        enhanced_vln_plot(Macrophages.All.subset, .x))

Macrophages.All.subset$module_group <- factor(Macrophages.All.subset$module_group ,
                                              levels = c("Module Low","Module High"))

single_vln_plot <- function(gene){

x <- VlnPlot(Macrophages.All.subset, 
                 group.by = "module_group", 
                 features = gene, cols = c("#cfcbd6","#6f45f6")) + 
  theme(text = element_text(size = 40)) +
  theme(axis.text.x = element_text(size = 40)) +
  theme(axis.text.y = element_text(size = 40)) + 
  theme(legend.position = "none") + 
  ylab("Normalized Expression") + 
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(face = "bold"))
  
  ggsave(x, file = glue("20201202_Expressors_Non_Expressors_ViolPlot_{gene}.png"), device = "png", 
         height = 20, width = 20)
  
} 

CCL3 <- VlnPlot(Macrophages.All.subset, 
                 group.by = "module_group", 
                 features = "CCL3", cols = c("#cfcbd6","#6f45f6")) + 
  theme(text = element_text(size = 40)) +
  theme(axis.text.x = element_text(size = 40)) +
  theme(axis.text.y = element_text(size = 40)) + 
  theme(legend.position = "none") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank())



x <- add_global_label((FABP4| CCL3), Xlab = "", size = 15, Xgap = .04)
ggsave(x, file = "20201015_Expressors_Non_Expressors_ViolPlot.png", device = "png", 
       height = 20, width = 20)


walk(c("CD68", "MRC1", "CCL4", "CXCL10", "CCL20", "CCL4L2"), single_vln_plot)


