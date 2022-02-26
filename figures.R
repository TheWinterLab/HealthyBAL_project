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


# Labeled Heatmap -----------------------------------------------

DimPlot(bal.all.control, label = TRUE, 
        reduction = "umap", 
        label.size = 7, 
        pt.size = 1,
        cols = c("#00a9ff","#cd9600","#00eea1","#e4ff00","#bf5e58","#c77cff","#ff61cc"))  +
  theme(legend.position = "none")

ggsave(filename = "20200721_Labeled_UMAP.pdf", 
       plot = last_plot(), 
       height = 10, 
       width = 10, 
       device = cairo_pdf)

#PNG

p1 <- DimPlot(bal.all.control, label = TRUE, reduction = "umap", label.size = 7,
              cols = c("#00a9ff","#cd9600","#00eea1","#e4ff00","#bf5e58","#c77cff","#ff61cc"))  +
  theme(legend.position = "none") # +
#ggtitle("Labeled Clusters")

CairoPNG("20200721_Labeled_UMAP.png", bg="transparent",
         height = 960, width = 960)
(
  p1
)
dev.off()

#Top 3 Genes Heatmap -----------------------------------------------------

bal.all.control.labeled.markers <- FindAllMarkers(object = bal.all.control, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(bal.all.control.labeled.markers %>% group_by(cluster), file="20190917_all_am_labeled_markers.txt", row.names=FALSE, sep="\t")
bal.all.control.labeled.markers <- read.delim(file="20190917_all_am_labeled_markers.txt")
cluster.averages <- AverageExpression(object = bal.all.control, return.seurat = TRUE)
top3 <- bal.all.control.labeled.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

CairoPDF("20201012_Top3_Genes_heatmap.pdf", bg="transparent",
         width = 8, height = 8)
(DoHeatmap(cluster.averages, features = top3$gene, 
           draw.lines = FALSE, size = 6, 
           group.colors = c("#00a9ff","#cd9600","#00eea1","#e4ff00","#bf5e58","#c77cff","#ff61cc"),
) +
    viridis::scale_fill_viridis() + 
    guides(fill = guide_legend(title="Relative Expression", reverse = TRUE, guides = FALSE),
           color = guide_colorbar())
)
dev.off()

# Umap with cell numbers --------------------------------------------------

umapCellAnno(bal.all.control, 
             #use.cols = c("#00a9ff","#cd9600","#00eea1","#e4ff00","#bf5e58","#c77cff","#ff61cc"),
             label.size = 7, point.size = .1, cell.legend.size = 12)

ggsave(filename = "20200825_Unlabeled_UMAP_with_Cell_Numbers.pdf",
       plot = last_plot(), 
       height = 10, 
       width = 11, 
       device = cairo_pdf)

# Clustering by Participant -----------------------------------------------
DimPlot(bal.all.control, group.by = "orig.ident") + 
  scale_color_discrete(name = "Individuals", 
                       labels = c("BAL01","BAL02", "BAL03","BAL04")) + 
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15)) + 
  guides(colour = guide_legend(override.aes = list(size = 3))) 

ggsave(filename = "20201014_Unlabeled_UMAP_Grouped_By_Cell_Identity.pdf",
       plot = last_plot(), 
       height = 10, 
       width = 11, 
       device = cairo_pdf)

# Cluster Correlation -----------------------------------------------------

cluster.averages <- AverageExpression(object = bal.all.control, show.progress = TRUE)
combined.corr <- cor(cluster.averages$SCT ,method="pearson")

ggsave(filename = "20200721_Cluster_Correlation.pdf",
       plot = pheatmap(combined.corr), 
       height = 6, 
       width = 6, 
       device = cairo_pdf)

# Unlabled UMAP with Cell Numbers -----------------------------------------
p1 <- umapCellAnno(bal.all.control, counts.as.title = TRUE, 
                   label.size = 17, point.size = .1,
                   axis.text.x.bottom.size  = 40,
                   axis.text.y.left.size = 40,
                   axis.title.y.size = 37,
                   axis.title.x.size = 37,
                   legend.text.size = 20)

CairoPDF(file = "20201014_Unlabeled_UMAP_with_Cell_Numbers.pdf",bg = "transparent",
         width = 20, height = 20)
(
  p1
)


dev.off()

CairoPNG(file = "20201014_Unlabeled_UMAP_with_Cell_Numbers.png",
         width = 1200, height = 1200)
(
  p1
)


dev.off()


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
         width = 10, height = 10)
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
         width = 10, height = 10)
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
         width = 10, height = 10)
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
         width = 10, height = 10)
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

CairoPDF("20200721_Unlabeled_UMAP.pdf", bg="transparent",
         width = 10, height = 10)
(
  p1
)
dev.off()





# Top 3 Gene Heatmap Unlabeled --------------------------------------------

bal.all.control.labeled.markers <- FindAllMarkers(object = bal.all.control, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(bal.all.control.labeled.markers %>% group_by(cluster), file="20190917_all_am_labeled_markers.txt", row.names=FALSE, sep="\t")
bal.all.control.labeled.markers <- read.delim(file="20190917_all_am_labeled_markers.txt")
cluster.averages <- AverageExpression(object = bal.all.control, return.seurat = TRUE)
top3 <- bal.all.control.labeled.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)

CairoPDF("20201014_Top3_Genes_Unlabeled_heatmap.pdf", bg="transparent",
         width = 20, height = 20)


(DoHeatmap(cluster.averages, features = top3$gene, size = 10,
           draw.lines = FALSE,
           group.colors = c("#00a9ff","#cd9600","#00eea1","#e4ff00","#bf5e58","#c77cff","#ff61cc")) + 
    viridis::scale_fill_viridis() + 
    theme(text = element_text(size = 25))+
    theme(axis.title.y = element_text(size = 30))+
    guides(fill = guide_legend(title="Relative Expression", reverse = TRUE, guides = FALSE),
           color = guide_colorbar())
)
dev.off()


# Unlabeled Cluster Correlation ---------------------------------------------------
library(circlize)
cluster.averages <- AverageExpression(object = bal.all.control, show.progress = TRUE)
combined.corr<-cor(cluster.averages$SCT,method="pearson")

heatmap.colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)

#pheatmap(combined.corr)

CairoPDF("20201111_Unlabeled_Cluster_Correlation.pdf", bg="transparent",
         width = 7, height = 6)
(
  Heatmap(combined.corr, col = heatmap.colors,
          show_column_names = FALSE, 
          heatmap_legend_param = list(title = "Correlation"),
          row_names_gp = gpar(fontface = "bold"),
          column_names_gp = gpar(fontface = "bold"),
          bottom_annotation = HeatmapAnnotation(
            text = anno_text(colnames(combined.corr), rot = 45, 
                             gp = gpar(fontface = "bold"))))
)

dev.off()


# Dot Plot ----------------------------------------------------------------

cell.table <- as_tibble(FetchData(bal.all.control, vars = c("orig.ident","Seurat_Assignment")))
count <- cell.table %>% group_by(orig.ident, Seurat_Assignment) %>% tally()
counts <- count %>% mutate(Percentage = case_when(
  orig.ident == "bal_ctrl_01" ~ n / (count %>% tally(wt = n))$n[[1]] * 100,
  orig.ident == "bal_ctrl_02" ~ n / (count %>% tally(wt = n))$n[[2]] * 100,
  orig.ident == "bal_ctrl_03" ~ n / (count %>% tally(wt = n))$n[[3]] * 100,
  orig.ident == "bal_ctrl_04" ~ n / (count %>% tally(wt = n))$n[[4]] * 100
))

p1 <- ggplot(counts, aes(x=Seurat_Assignment, y=Percentage)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize = .5) + 
  scale_y_continuous(trans='log2', breaks = c(.1,.5,1,2,5,10,20,50,80,100),limits = c(.1,100)) +
  labs(x = "Cell Type") + 
  scale_fill_discrete(name = "Individual", labels = c("BAL01", "BAL02", "BAL03","BAL04")) + 
  theme(text = element_text(size = 50))  + 
  geom_boxplot(alpha = .1) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

CairoPDF("20200728_CellType_Percentage_by_Individual.pdf", bg="transparent",
         width = 20, height = 20)
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
                     labels = c(.050,.100,.500,1.00,2.00,5.00,10.0,20.0,50.0,80.0,100),
                     guide = guide_axis(check.overlap = TRUE))+ 
  theme(axis.title = element_text(size = 45),
        legend.text = element_text(size = 45),
        legend.title = element_text(size = 45),
        axis.text.x = element_text(size = 45, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 45),
        panel.background = element_rect(fill = "#FFFFFF"),
        panel.border = element_rect(fill = NA)) + 
  labs(x = "Cell Type",
       y = "% of Cells Annotated by Cell Type") +
  theme(legend.key.size = unit(50,"point"))



CairoPDF("20201012_CellType_Percentage_by_Individual_With_Color.pdf", bg="transparent",
         width = 20, height = 20)
(
  p2
)

dev.off()



# Volcano Plot ------------------------------------------------------------

load("20200729_bal_mp.Robj")
cluster.names = c("Macrophages","MIP-1 Macrophages")
names(cluster.names) <- levels(x = bal.mp)
bal.mp <- Seurat::RenameIdents(bal.mp, cluster.names)
bal.mp@meta.data$Seurat_Assignment <- Idents(bal.mp)


bal.mp <- FindVariableFeatures(bal.mp, nfeatures = 50)
p1 <- VariableFeaturePlot(bal.mp, pt.size = 2)

p2 <- LabelPoints(p1, points= head(VariableFeatures(bal.mp), 10), size = 11, repel = T) + 
  theme(text = element_text(size = 40),
        axis.ticks = element_line(size = 3),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        panel.border = element_rect(fill = "NA"),
        legend.position = "none")

CairoPDF("20201014_Macrophage_Variable_Gene_Plot.pdf", bg="transparent",
         width = 20, height = 20)
(
  p2
)

dev.off()



# Volcano Plot Separated by Participant ----------------------------------------------------
load("20200729_bal_mp.Robj")
library(magrittr); library(glue);
cluster.names = c("Macrophages","MIP-1 Macrophages")
names(cluster.names) <- levels(x = bal.mp)
bal.mp <- Seurat::RenameIdents(bal.mp, cluster.names)
bal.mp@meta.data$Seurat_Assignment <- Idents(bal.mp)

split.bal.obj <- SplitObject(bal.mp, split.by = "orig.ident")
split.bal.obj %<>% map(~ FindVariableFeatures(.x, nfeatures = 50))
p1 <- map(split.bal.obj, ~ VariableFeaturePlot(.x, pt.size = 2))

mip1_count <- map(split.bal.obj, ~ FetchData(.x, "ident") %>% 
                    as_tibble(rownames = "ID") %>%
                    group_by(ident) %>% 
                    tally() %>% 
                    pull(n) %>% 
                    .[2])

plot_volc <- function(p1, bal_mp) { 
  LabelPoints(p1, points = head(VariableFeatures(bal_mp), 10), 
              size = 8, repel = T) + 
    theme(text = element_text(size = 20),
          axis.ticks = element_line(size = 3),
          axis.title.x = element_text(size = 30),
          axis.title.y = element_text(size = 30),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          panel.border = element_rect(fill = "NA"),
          legend.position = "none")
} 

x <- map2(p1, split.bal.obj, plot_volc)

# x[[1]] <- x[[1]] + annotate(geom = "text", x = 200, y = 20, 
#                             label = glue("MIP-1 Macrophages: \n {mip1_count[[1]]}"),
#                             size = 10,
#                             color = "red") 
# 
# x[[2]] <- x[[2]] + annotate(geom = "text", x = 200, y = 20, 
#                             label = glue("MIP-1 Macrophages: \n {mip1_count[[2]]}"),
#                             size = 10,
#                             color = "red") 
# 
# x[[3]] <- x[[3]] + annotate(geom = "text", x = 200, y = 20, 
#                             label = glue("MIP-1 Macrophages: \n {mip1_count[[3]]}"),
#                             size = 10,
#                             color = "red") 
# 
# x[[4]] <- x[[4]] + annotate(geom = "text", x = 200, y = 20, 
#                             label = glue("MIP-1 Macrophages: \n {mip1_count[[4]]}"),
#                             size = 10,
#                             color = "red") 


file.names <- paste0("20201202_Macrophage_Variable_Gene_Plot_", c("BAL01.pdf","BAL02.pdf",
                                                                  "BAL03.pdf","BAL04.pdf"))
for (i in seq_along(x)) { 
  
  CairoPDF(file.names[[i]], bg="transparent",
           width = 20, height = 20)
  
  (
    print(x[[i]])
  )
  dev.off()
}

file.names <- paste0("20201202_Macrophage_Variable_Gene_Plot_", c("BAL01.png","BAL02.png",
                                                                  "BAL03.png","BAL04.png"))

for (i in seq_along(x)) { 
  
  CairoPNG(file.names[[i]], bg="transparent",
           width = 1200, height = 1200)
  
  (
    print(x[[i]])
  )
  dev.off()
}



# Bal MP Grouped by Participant Stacked Bar Plot --------------------------------------------
load("20200729_bal_mp.Robj")
library(magrittr); library(glue);
cluster.names = c("Macrophages","MIP-1 Macrophages")
names(cluster.names) <- levels(x = bal.mp)
bal.mp <- Seurat::RenameIdents(bal.mp, cluster.names)
bal.mp@meta.data$Seurat_Assignment <- Idents(bal.mp)

mip1_count <- FetchData(bal.mp, c("orig.ident", "Seurat_Assignment")) %>% as_tibble(rownames = "ID") %>%
                    group_by(`orig.ident`, Seurat_Assignment) %>% tally()



mip_sbp <- ggplot(mip1_count, aes(orig.ident, n, fill = Seurat_Assignment )) + 
  geom_bar(position = "stack", stat = "identity") +
  theme_bw() + 
  theme(legend.title = element_blank(),
        text = element_text(size = 40),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "#000000")) + 
  scale_x_discrete(labels = c("BAL01", "BAL02", "BAL03", "BAL04")) +
  xlab("Individual") + 
  ylab("Cells") 

mip_sbp <- mip_sbp + annotate(geom = "text", x = c("bal_ctrl_01", "bal_ctrl_02", "bal_ctrl_03", "bal_ctrl_04"),
                             y = c(8000, 9000, 6500, 8000), label = mip1_count$n[c(1,3,5,7)], size = 12)


# mip_sbp <- mip_sbp + annotate(geom = "text", x = c("bal_ctrl_01", "bal_ctrl_02", "bal_ctrl_03", "bal_ctrl_04"),
#                    y = c(300, 400, 200, 200), label = mip1_count$n[c(2,4,6,8)], size = 10)

mip_sbp <- mip_sbp + annotate(geom = "text", x = c("bal_ctrl_01", "bal_ctrl_02", "bal_ctrl_03", "bal_ctrl_04"),
                              y = c(200, 200, 200, 200), label = mip1_count$n[c(2,4,6,8)], size = 12)


CairoPDF("20201202_MIP1_Stacked_Barplot.pdf", bg="transparent",
         width = 25, height = 20)

( 
  mip_sbp
) 


dev.off()


CairoPNG("20201202_MIP1_Stacked_Barplot.png", bg="transparent",
         width = 1600, height = 1200)


( 
  mip_sbp
) 


dev.off()

# Feature Plot Split by Participant Macrophages Only ---------------------------------------
load("20200729_bal_mp.Robj")
library(magrittr); library(glue);
cluster.names = c("Macrophages","MIP-1 Macrophages")
names(cluster.names) <- levels(x = bal.mp)
bal.mp <- Seurat::RenameIdents(bal.mp, cluster.names)
bal.mp@meta.data$Seurat_Assignment <- Idents(bal.mp)

split.bal.obj <- SplitObject(bal.mp, split.by = "orig.ident")
features <- c("CD68", "MRC1", "CCL4", "CXCL10", "CCL20", "CCL4L2") 

fp <- function(feature){ 
  x <- map(split.bal.obj, ~ FeaturePlot(.x, features = feature),
           ncol = 1, cols = c("grey", "red"), pt.size = .1)
  
  file.name <- paste0(feature, ".png")
  file.name_2 <- map_chr(c("BAL01","BAL02","BAL03","BAL04"), ~ 
                           glue::glue("20201118_FeaturePlot_{.x}"))
  file.name <- paste(file.name_2, file.name, sep = "_")
  
  for (i in seq_along(x)){ 
    CairoPNG(file.name[[i]], bg="transparent",
             width = 375, height = 375)
    (
      print(x[[i]])
    )
    dev.off()
  }
} 

# MP Feature Plot CD68, MRC1, CXCL10, ... --------------------------------------------------
features <- c("CD68", "MRC1", "CCL4", "CXCL10", "CCL20", "CCL4L2") 

fp_vector <- vector(mode = "list", length = length(features))
for (i in seq_along(fp_vector)){
  
fp_vector[[i]] <- FeaturePlot(bal.mp, features = features[[i]],
    ncol = 1, cols = c("grey", "red"), pt.size = .1)
}

file.name <- paste0(features, ".png")
file.name <- paste("20201202_BAL_MP_FeaturePlot", file.name, sep = "_")

for (i in seq_along(fp_mp)) {
  CairoPNG(file.name[[i]], bg="transparent",
         width = 375, height = 375)
  (
    print(fp_vector[[i]])
  )
  dev.off()
}

# MP Module Score ---------------------------------------------------------


bal.mp <- AddModuleScore(object = bal.mp,features = list(c("CCL3", "CCL4", "CXCL10", "CCL4L2" ,"CCL20")),
                         name = "Module_Scores")

CairoPDF("20200825_Macrophage_MIP_Subcluster_Module_Score.pdf", bg="transparent",
         width = 20, height = 20)

(
  FeaturePlot(bal.mp, features = "Module_Scores1", pt.size = 2)
  + ggtitle("MIP-1 Module Score") + 
    theme(text = element_text(size = 20)) + 
    theme(axis.text = element_text(size = 20))
)

dev.off()



# Variable Gene Matrix ----------------------------------------------------

cor.mat <- variableGeneMatrix(bal.mp, variable.genes = 50, return.table = TRUE)
cor.mat.plot <- pheatmap(cor.mat, angle_col = 90, fontsize = 13)


CairoPDF("20200729_Variable_Gene_Correlation_Matrix.pdf", bg="transparent",
         width = 20, height = 20)
(
  cor.mat.plot
)

dev.off()




# Macrophage UMAP ---------------------------------------------------------

p1 <- umapCellAnno(bal.mp)

CairoPDF(file = "20200728_Labeled_Macrophages_UMAP_with_Cell_Numbers.pdf",bg = "transparent",
         width = 11, height = 10)


(
  umapCellAnno(bal.mp, 
               label.size = 7, point.size = .1, cell.legend.size = 12)
)

dev.off()




# Differential Genes MIP and Macrophages ----------------------------------
# Differential Genes MIP and Macrophages/Variable Gene Heatmap ------------

mip.all.markers <- FindMarkers(bal.mp, ident.1 = "MIP-1 Macrophages", ident.2 = "Macrophages")
mip.markers <- cbind(gene = row.names(mip.all.markers), mip.all.markers)

top5.mip.markers <- mip.markers %>% top_n(n = 10, wt = avg_logFC)
bot5.mip.markers <- mip.markers %>% top_n(n = 10, wt = -avg_logFC)
am.cells <- colnames(x = subset(x = bal.mp,
                                idents = c("MIP-1 Macrophages", "Macrophages"),
                                downsample = 500))


CairoPDF("20201014_Macrophage_MIP_Var_Gene_Heatmap_tophalf.pdf", bg="transparent",
         width = 20, height = 20)

(
  DoHeatmap(bal.mp, 
            cells = am.cells, 
            features = rownames(top5.mip.markers), 
            size = 12,
            angle = 0,
            hjust = .5) +
    scale_fill_viridis() +
    theme(axis.text.y = element_text(size = 25),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20)) +  
    guides(fill = guide_legend(title="Relative Expression", reverse = TRUE, guides = FALSE),
           color = guide_colorbar())
) 



dev.off()

CairoPDF("20200914_Macrophage_MIP_Var_Gene_Heatmap_bothalf.pdf", bg="transparent",
         width = 20, height = 20)

(
  DoHeatmap(bal.mp, 
            cells = am.cells, 
            features = rownames(bot5.mip.markers), 
            angle = 0,
            hjust = .5,
            size = 10) + 
    scale_fill_viridis() +
    theme(text = element_text(size = 20))  +
    guides(colour = FALSE) 
) 

dev.off()




# Over Representation Bar Plot --------------------------------------------


all.genes <- read_tsv("20191218_mip1_gsea.rnk", col_names = c("gene_name","fold_change"))
mip.markers <- mip.marker$gene
entrez.all.genes <- mapIds(org.Hs.eg.db, all.genes$gene_name, 'ENTREZID', 'SYMBOL')
entrez.mip.markers <- mapIds(org.Hs.eg.db, mip.markers, 'ENTREZID', 'SYMBOL')


ego2 <- enrichGO(gene = entrez.mip.markers,
                 OrgDb  = org.Hs.eg.db,
                 universe = entrez.all.genes,
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

barplot(ego2, showCategory = 20)


CairoPDF("20200729_Barplot_Over_Representation_Test.pdf", bg="transparent",
         width = 10, height = 10)

(barplot(ego2, showCategory = 20))

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

CairoPDF("202000806_CellType_Percentage_BAL05_With_Color.pdf", bg="transparent",
         width = 20, height = 20)
(
  p2
)

dev.off()



# Ambient RNA Scatter Plot ------------------------------------------------

callout_genes <- c("FABP4", "CCL5", "LTB", "MS4A1", "CTSB", "SCGB3A1", "CAPS", "TOP2A")
sample.names <- c("BAL01","BAL02","BAL03","BAL04")
sample.scatter <- map(sample.names, sampleScatter)
sample.scatter2 <- map2(sample.scatter, sample.names, ~ plot.grid(.x,.y))

CairoPDF("202000807_Ambient_RNA_Scatter_Plot.pdf", bg="transparent",
         width = 20, height = 20)
(
  cowplot::plot_grid(sample.scatter2[[1]],sample.scatter2[[2]], sample.scatter2[[3]], sample.scatter2[[4]],
                     nrow = 2, ncol =2)
) 

dev.off()


matSamp <- extractMatrixSums(c("BAL01","BAL02","BAL03","BAL04"))

celltypeGraph(unique(bal.all.control$orig.ident), matSamp, bal.all.control)



# HL Donor Macrophages Module Scores --------------------------------------
load("20200811_Macrophages_Donors_Only.Robj")
Macrophages.All.subset <- RunUMAP(Macrophages.All.subset, dims = 1:30, seed.use = 24)

Macrophages.All.subset <- AddModuleScore(object = Macrophages.All.subset,
                                         features = list(c("CCL3", "CCL4", "CXCL10", "CCL4L2" ,"CCL20")),
                                         name = "Module_Scores")


CairoPDF("20200825_HL_Macrophage_Donor_Module_Score.pdf", bg="transparent",
         width = 20, height = 20)

(
  FeaturePlot(Macrophages.All.subset, features = "Module_Scores1", pt.size = 2,
              label = TRUE, label.size = 20) + ggtitle("MIP-1 Module Score") + 
    theme(text = element_text(size = 20)) + 
    theme(axis.text = element_text(size = 20))
)

dev.off()


# HL Donor Macrophages Feature Plot ---------------------------------------

CairoPDF("20200818_HL_Macrophages_FeaturePlot_FABP4.pdf", bg="transparent",
         width = 20, height = 20)

(FeaturePlot(Macrophages.All.subset, features = c("FABP4"),
             ncol = 1, cols = c("grey", "red"), coord.fixed = TRUE, pt.size = .1)
)

dev.off()


# HL UMAP -----------------------------------------------------------------


p1 <- DimPlot(Macrophages.All.subset, label = TRUE, reduction = "umap", label.size = 12)  +
  theme(legend.position = "none") 

CairoPDF("20200826_Unlabeled_UMAP.pdf", bg="transparent",
         width = 10, height = 10)
(
  p1
)
dev.off()


# Seurat Dot Plot ---------------------------------------------------------

#Read in clusters to help manually change y-axis values:
cluster <- read_csv("cluster_csv.csv")

CairoPDF("20200916_Seurat_DotPlot_Celltype_Grouped_by_Individual.pdf", bg="transparent",
         width = 20, height = 20)

p1 <- DotPlot(bal.all.control, 
              features = c("FABP4","CCL5","MS4A1", "CTSB","SCGB3A1","CAPS","TOP2A"),
              group.by = "Seurat_Assignment",
              split.by = "orig.ident",
              dot.scale = 15, cols = c("#f42013","#85b312","#08aeb1","#910df2")) +
  theme(panel.background = element_rect(fill = "NA"),
        text = element_text(size = 25),
        axis.title = element_blank()) +
  theme(panel.border = element_rect(color = "#000000", fill = "NA")) +
  scale_y_discrete(labels = rep(4:1, 8)) +
  coord_fixed(ratio = 1/2)


features = c("FABP4","CCL5","MS4A1", "CTSB","SCGB3A1","CAPS","TOP2A")

x <- p1$data %>%
  as_tibble() %>%
  unite(col = identifier, c("id", "features.plot"), sep = "_", remove = FALSE) %>%
  dplyr::select(-id) %>%
  mutate(cont_identifier = 1:196) %>%
  dplyr::rename(`Percent Expressed` = pct.exp) 

cells <- c("Dividing Cells_", 
           "Ciliated Cells_", 
           "Club Cells_", 
           "Dendritic Cells_", 
           "B Cells_", 
           "T/NK Cells_", 
           "Macrophages_")

cells <- rep(cells, each = 28)
bctrl <- rep(rep(paste0("bal_ctrl_0", 1:4), each = 7), 7)
cells <- paste0(cells, bctrl)
cells <- paste(cells, unique(x$features.plot), sep = "_")


x %<>% arrange(factor(identifier, levels = cells)) %>%
  mutate(y = rep(1:28, each = 7))


y_axis <- str_c(rep(4:1, 7)   , "\n")
y_axis[y_axis == "4\n"] <- "4\n---"

new <- ggplot(data = x, aes(x = features.plot, 
                            y = y,
                            size = `Percent Expressed`)) +
  geom_point(color = x$colors) + 
  scale_size_continuous(range = c(minSize = 1, maxSize = 10)) +
  theme(panel.background = element_rect(fill = "NA"),
        text = element_text(size = 12),
        axis.title = element_blank(),
        plot.margin = margin(2,2,2,4,"cm"),
        axis.ticks.y.right = element_blank(),
        axis.text.y.left = element_text(size = 13)) + 
  #axis.text.y.right = element_text(angle = 90)) +
  theme(panel.border = element_rect(color = "#000000", fill = "NA")) + 
  scale_y_continuous(breaks = 1:28, labels = rep(4:1,7),
                     sec.axis = sec_axis(~., 
                                         breaks = c(.8, 4.8, 8.8, 12.8, 16.8, 20.8,  24.8), 
                                         labels = c("---------------",
                                                    "---------------",
                                                    "---------------",
                                                    "---------------",
                                                    "---------------",
                                                    "---------------",
                                                    "---------------"))) + 
  labs(tag = "BAL #") + theme(plot.tag.position = c(-.01, 1),
                              plot.tag = element_text(face = "plain", size = 12)) +
  scale_size_continuous(name = "")

web_browse(new, width = 9, height = 11)


CairoPDF("20201112_Seurat_DotPlot_Celltype_Grouped_by_Individual.pdf", bg="transparent",
         width = 9, height = 10)

(
  new
)

dev.off()




(
  DotPlot(bal.all.control, 
          features = c("FABP4","CCL5","MS4A1", "CTSB","SCGB3A1","CAPS","TOP2A"),
          group.by = "Seurat_Assignment",
          split.by = "orig.ident",
          dot.scale = 15, cols = c("#f42013","#85b312","#08aeb1","#910df2")) +
    theme(panel.background = element_rect(fill = "NA")) + 
    theme(text = element_text(size = 25)) +
    scale_y_discrete(labels = rev(cluster)) +
    theme(panel.border = element_rect(color = "#000000", fill = "NA")) +
    ylab("Cell Types - Grouped by Individual") + 
    xlab("Genes") 
) 

dev.off()

# Ambient RNA Excluded/Included Four by Four -----------------------------------------------

sample.matrix <- extractMatrixSums(c("BAL01","BAL02","BAL03","BAL04"))
mat <- fourByFour(c("bal_ctrl_01","bal_ctrl_02","bal_ctrl_03","bal_ctrl_04"), 
                  sample.matrix = sample.matrix)

dimnames(mat) <- list(c("BAL01","BAL02","BAL03","BAL4"), c("BAL01","BAL02","BAL03","BAL04"))

heatmap.colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)


CairoPDF("20210126_Ambient_RNA_Excluded_Included_By_Individual_Heatmap.pdf", bg="transparent",
         width = 6, height = 6)

(

  Heatmap(mat, col = heatmap.colors,
               cluster_rows = FALSE,
          show_row_names = TRUE,
          show_column_names = TRUE,
               cluster_columns = FALSE,
               column_title = "Included",
               row_title = "Excluded",
          heatmap_legend_param = list(title = ""))
  
)


dev.off()


#(
#  pheatmap(mat, fontsize = 10, cluster_cols = FALSE, cluster_rows = FALSE, border_color = NA)
#)

dev.off()

heatmap.colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)


CairoPDF("20201112_Unlabeled_Cluster_Correlation.pdf", bg="transparent",
         width = 7, height = 6)
(
  x <- Heatmap(mat, col = heatmap.colors,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               heatmap_legend_param = list(title = "Correlation"))
  
)

dev.off()

