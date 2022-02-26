library(Seurat)
library(pheatmap)
library(ComplexHeatmap)
library(Matrix)
library(Cairo)
library(tidyverse)
library(data.table)
library(SingleR)
library(ggpubr)
library(EnhancedVolcano)
library(ggrepel)
library(viridis)
library(org.Hs.eg.db)
library(RColorBrewer)
library(glue)
library(magrittr)
library(clusterProfiler)
library(data.table)

# Functions that end with a 2 are enhanced versions of their predecessor


renameClusters_2 <- function(seurat.obj, 
                           cluster.names,
                           annotation.name = NULL) { 
  
  if(is.null(annotation.name)){ 
    annotation.name <- "Seurat_Assignment"
  }
  
  names(cluster.names) <- levels(x = seurat.obj)
  seurat.obj <- RenameIdents(seurat.obj, cluster.names)
  seurat.obj@meta.data[[annotation.name]] <- Idents(seurat.obj)
  return(seurat.obj)
  
}



extractMatrixSums <- function(samples) { 
  matrix_list_raw <- vector(mode = "list", length = length(samples))
  for (i in 1:length(samples)) { 
    barcode_list_raw <- vector(mode = "list", length = length(samples))
    barcode_list_filtered <- vector(mode = "list", length = length(samples))
    features <- read.delim("features.tsv.gz",header = FALSE, stringsAsFactors = FALSE)
    matrix_list_raw[[i]] <- readMM(paste(samples[[i]],"_raw_matrix.mtx.gz",sep=""))
    barcode_list_raw[[i]] <- read.delim(paste(samples[[i]],"_raw_barcodes.tsv.gz",sep=""),header = FALSE, stringsAsFactors = FALSE)
    barcode_list_filtered[[i]] <- read.delim(paste(samples[[i]],"_filtered_barcodes.tsv.gz",sep=""),header = FALSE, stringsAsFactors = FALSE)
    colnames(matrix_list_raw[[i]]) <- barcode_list_raw[[i]]$V1
    rownames(matrix_list_raw[[i]]) <- features$V2
    barcode_list_raw[[i]] <- setdiff(barcode_list_raw[[i]]$V1, barcode_list_filtered[[i]]$V1)
    matrix_list_raw[[i]] <- matrix_list_raw[[i]][,barcode_list_raw[[i]]]
    matrix_list_raw[[i]] <- as.data.frame(as.table(rowSums(matrix_list_raw[[i]])))
    colnames(matrix_list_raw[[i]]) <- c("gene_name","counts")
  }
  return(matrix_list_raw)
}


extractMatrixSums2 <- function(sample.name) { 
  features <- read_tsv("features.tsv.gz", col_names = FALSE, col_types = cols())$X2
  raw_barcode <- read_tsv(paste0(sample.name,"_raw_barcodes.tsv.gz"), col_names = FALSE, col_types = cols())$X1
  filtered_barcode <- read_tsv(paste0(sample.name,"_filtered_barcodes.tsv.gz"), col_names = FALSE, col_types = cols())$X1
  raw_counts <- readMM(paste(sample.name,"_raw_matrix.mtx.gz",sep=""))
  colnames(raw_counts) <- raw_barcode
  raw_counts <- raw_counts[,setdiff(raw_barcode, filtered_barcode)]
  raw_counts <- tibble(gene_name = features, counts = rowSums(raw_counts))
  return(raw_counts)
}



sampleScatter <- function(i) { 
  features <- read_tsv("features.tsv.gz", col_names = FALSE)$X2
  matrix_list_raw <- readMM(paste(i,"_raw_matrix.mtx.gz",sep=""))
  matrix_list_filtered <- readMM(paste(i,"_filtered_matrix.mtx.gz",sep=""))
  barcode_list_raw <- read_tsv(paste(i,"_raw_barcodes.tsv.gz",sep=""), col_names = FALSE)$X1
  barcode_list_filtered<- read_tsv(paste(i,"_filtered_barcodes.tsv.gz",sep=""), col_names = FALSE)$X1
  colnames(matrix_list_raw) <- barcode_list_raw
  rownames(matrix_list_raw) <- features
  colnames(matrix_list_filtered) <- barcode_list_filtered
  rownames(matrix_list_filtered) = features
  barcode_list_raw <- setdiff(barcode_list_raw, barcode_list_filtered)
  matrix_list_raw <- matrix_list_raw[,barcode_list_raw]
  mat_combined <- as.data.frame(cbind(rowSums(matrix_list_raw),rowSums(matrix_list_filtered)))
  setDT(mat_combined, keep.rownames = TRUE)
  colnames(mat_combined) <- c("gene_name", "excluded_total_counts","included_total_counts")
  excluded_CPM <- sum(mat_combined$excluded_total_counts)
  included_CPM <- sum(mat_combined$included_total_counts)
  mat_combined$excluded_total_counts <- (mat_combined$excluded_total_counts/excluded_CPM)*10^6
  mat_combined$included_total_counts <- (mat_combined$included_total_counts/included_CPM)*10^6
  #excluded_counts <- mat_combined[(mat_combined$included_total_counts>0 & mat_combined$excluded_total_counts>0),]
  return(mat_combined)
} 

excluded_included_heatmap <- function(i) { 
  features <- read_tsv("features.tsv.gz", col_names = FALSE)$X2
  matrix_list_raw <- readMM(paste(i,"_raw_matrix.mtx.gz",sep=""))
  matrix_list_filtered <- readMM(paste(i,"_filtered_matrix.mtx.gz",sep=""))
  barcode_list_raw <- read_tsv(paste(i,"_raw_barcodes.tsv.gz",sep=""), col_names = FALSE)$X1
  barcode_list_filtered<- read_tsv(paste(i,"_filtered_barcodes.tsv.gz",sep=""), col_names = FALSE)$X1
  colnames(matrix_list_raw) <- barcode_list_raw
  rownames(matrix_list_raw) <- features
  colnames(matrix_list_filtered) <- barcode_list_filtered
  rownames(matrix_list_filtered) = features
  barcode_list_raw <- setdiff(barcode_list_raw, barcode_list_filtered)
  matrix_list_raw <- matrix_list_raw[,barcode_list_raw]
  mat_combined <- as.data.frame(cbind(rowSums(matrix_list_raw),rowSums(matrix_list_filtered)))
  setDT(mat_combined, keep.rownames = TRUE)
  colnames(mat_combined) <- c("gene_name", "excluded_total_counts","included_total_counts")
  excluded_CPM <- sum(mat_combined$excluded_total_counts)
  included_CPM <- sum(mat_combined$included_total_counts)
  mat_combined$excluded_total_counts <- (mat_combined$excluded_total_counts/excluded_CPM)*10^6
  mat_combined$included_total_counts <- (mat_combined$included_total_counts/included_CPM)*10^6
  colnames(mat_combined) <- c("gene_name", paste(i,"excluded_total_counts",sep = "_"), paste(i,"included_total_counts",sep ="_"))
  return(mat_combined)
  
} 

plot.grid <- function(x, title) {
  ggplot(x,aes(x = included_total_counts, y = excluded_total_counts)) +
    geom_point(size = 1) +
    scale_x_continuous(trans = 'log2') +
    scale_y_continuous(trans = 'log2') +
    labs(x = "Included Total counts", y = "Excluded Total Counts") +
    geom_label_repel(size = 10, data = subset(x, gene_name %in% callout_genes),
                     aes(x = included_total_counts, y = excluded_total_counts, label = gene_name),
                     box.padding   = 0.35,
                     point.padding = 0.5,
                     segment.color = 'red',
                     segment.size = 1,
                     label.size = .75) +
    stat_cor(size = 10, method = "pearson") +
    theme(text = element_text(size = 28),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          panel.border = element_rect(fill = "NA"),
          panel.background = element_rect(fill = "#ffffff")) + 
    ggtitle(title) 
  
} 


hier <- function(annotation.name = "Seurat_Assignment",
                 down.sample = 50,
                 variable.genes = 50,
                 cluster,
                 seed = 1){
  
  set.seed(seed)
  cell.extract <- as_tibble(FetchData(bal.mp, vars = annotation.name),rownames = "Cells") %>%
    filter(!!as.name(annotation.name) == cluster) %>%
    slice_sample(n = down.sample, replace = FALSE)
  counts <- as_tibble(GetAssayData(bal.mp)[VariableFeatures(bal.mp)[seq(variable.genes)], cell.extract$Cells], rownames = "Genes") %>%
    pivot_longer(cols = -(Genes), names_to = "Cells") %>%
    pivot_wider(names_from = Genes, values_from = value) %>%
    right_join(y = cell.extract, by = "Cells")
  mat <- as.matrix(counts %>% select(-c(Cells,Seurat_Assignment)))
  rownames(mat) <- counts$Cells
  colnames(mat) <- colnames(counts %>% select(-c(Cells,Seurat_Assignment)))
  dist <- hclust(dist(mat))
  fviz_dend(dist)
  return(list(mat, dist))
} 


celltypeMatrix <- function(sample.list,
                           matrix.sample.list, 
                           seurat) { 
  clusters <- unique(seurat@meta.data$Seurat_Assignment)
  matrix.sum <- extractMatrixSums(matrix.sample.list)
  empty.matrix <- matrix(nrow = length(clusters) * length(sample.list), ncol = length(sample.list))
  sample.subset.merged <- c()
  row <- 1
  for (i in 1:length(clusters)) { 
    for (j in 1:length(sample.list)) { 
      for (k in 1:length(sample.list)) { 
        sample.subset <- setDT(FetchData(seurat, vars = c("orig.ident","Seurat_Assignment")), keep.rownames = TRUE)
        sample.subset <- sample.subset[sample.subset$orig.ident == sample.list[[k]] & sample.subset$Seurat_Assignment == clusters[[i]]]
        sample.subset <- setDT(as.data.frame(rowSums(GetAssayData(subset(seurat, cells = sample.subset$rn), assay = "RNA"))), keep.rownames = TRUE)
        colnames(sample.subset) <- c("gene_name", "counts")
        sample.subset.merged <- full_join(matrix.sum[[j]], sample.subset, by = "gene_name")
        sample.subset.merged[is.na(sample.subset.merged)] = 0
        sample.subset.merged$counts.x <- (sample.subset.merged$counts.x/sum(sample.subset.merged$counts.x))*10^6
        sample.subset.merged$counts.y <- (sample.subset.merged$counts.y/sum(sample.subset.merged$counts.y))*10^6
        sample.subset.merged <- sample.subset.merged[sample.subset.merged$counts.x > 0 & sample.subset.merged$counts.y >0,]
        sample.subset.merged <- log(sample.subset.merged[-1])
        empty.matrix[row,k] <- cor(sample.subset.merged$counts.x,sample.subset.merged$counts.y)
      }
      row <- row + 1
    }
  }
  return(empty.matrix)
}


celltypeGraph <- function(sample.list, matrix.sample.list, seurat) { 
  clusters <- unique(seurat@meta.data$Seurat_Assignment)
  clusters_print <- str_replace_all(clusters, "/","_")
  clusters_print <- str_replace_all(clusters_print, " ","_")
  sample.subset <- vector(mode = "list", length = length(matrix.sample.list) * length(clusters))
  for (i in 1:length(matrix.sample.list)) { 
    for (j in 1:length(clusters)) { 
      sample.subset <- setDT(FetchData(seurat, vars = c("orig.ident","Seurat_Assignment")), keep.rownames = TRUE)
      sample.subset <- sample.subset[sample.subset$orig.ident == sample.list[[i]]]
      sample.subset <- sample.subset[sample.subset$Seurat_Assignment == clusters[[j]]]
      sample.subset <- setDT(as.data.frame(rowSums(GetAssayData(subset(seurat, cells = sample.subset$rn), assay = "RNA"))), keep.rownames = TRUE)
      colnames(sample.subset) <- c("gene_name", "counts")
      sample.subset <- full_join(matrix.sample.list[[i]], sample.subset, by = "gene_name")
      sample.subset[is.na(sample.subset)] = 0
      colnames(sample.subset) <- c("gene_name","excluded_total_counts","included_total_counts")
      excluded_counts <- sample.subset[(sample.subset$included_total_counts>0 & sample.subset$excluded_total_counts>0),]
      excluded_CPM <- sum(excluded_counts$excluded_total_counts)
      included_CPM <- sum(excluded_counts$included_total_counts)
      excluded_counts$excluded_total_counts <- (excluded_counts$excluded_total_counts/excluded_CPM)*10^6
      excluded_counts$included_total_counts <- (excluded_counts$included_total_counts/included_CPM)*10^6
      callout_genes <- c("FABP4", "CCL5", "LTB", "MS4A1", "CTSB", "SCGB3A1", "CAPS", "TOP2A")
      ggplot_1 <- ggplot(excluded_counts, aes(x=included_total_counts, y=excluded_total_counts)) +
        geom_point(size = .375) +
        scale_x_continuous(trans='log2') +
        scale_y_continuous(trans='log2') +
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              panel.border = element_rect(fill = "NA"),
              panel.background = element_rect(fill = "#ffffff")) + 
        labs(x = "Included Total counts", y = "Excluded Total Counts") +
        geom_label_repel(data=subset(excluded_counts, gene_name %in% callout_genes),
                         aes(x=included_total_counts, y=excluded_total_counts, label = gene_name),
                         size = 4.75,
                         box.padding   = 0.35,
                         segment.size = .25,
                         point.padding = 0.2,
                         label.size = .70,
                         segment.color = 'red')
      ggplot_1 = ggplot_1 + stat_cor(size = 5.75, method = "pearson") 
      
     # s.l <- str_to_upper(gsub(pattern = "_ctrl_", replacement = "",sample.list[[i]]))
     # ggplot_1 = ggplot_1+ ggtitle(paste(s.l, clusters[[j]], sep = " ")) 

     ggsave(paste(paste(sample.list[[i]], clusters_print[j], sep = "_"),"_scatter_plot.png",sep=""), plot = ggplot_1 , width = 5, height = 5, dpi = 90, device = "png")
      
    }
    
    
  }
}


fourByFour <- function(original.ident, sample.matrix)
{ 
  empty.matrix <- matrix(nrow = length(original.ident), ncol = length(sample.matrix))
  sample.subset <- vector(mode = "list", length = length(sample.matrix))
  for (i in 1:length(original.ident))
  { 
    sample.subset[[i]] <- setDT(FetchData(bal.all.control, vars = "orig.ident"), keep.rownames = TRUE)
    sample.subset[[i]] <- sample.subset[[i]][sample.subset[[i]]$orig.ident == original.ident[[i]]]
    sample.subset[[i]] <- setDT(as.data.frame(rowSums(GetAssayData(subset(bal.all.control, cells = sample.subset[[i]]$rn), assay = "RNA"))), keep.rownames = TRUE)
    colnames(sample.subset[[i]]) <- c("gene_name", "counts")
  }
  
  for (i in 1:length(sample.matrix))
  {
    for (j in 1:length(sample.subset))
    {
      sample.subset.merged <- c()
      sample.subset.merged <- full_join(sample.matrix[[j]], sample.subset[[i]], by = "gene_name")
      sample.subset.merged[is.na(sample.subset.merged)] = 0
      sample.subset.merged <- sample.subset.merged[sample.subset.merged$counts.x > 0 & sample.subset.merged$counts.y >0,]
      sample.subset.merged <- log(sample.subset.merged[-1])
      empty.matrix[i,j] <- cor(sample.subset.merged$counts.x,sample.subset.merged$counts.y)
      
    }
  }
  return(empty.matrix)
}

extractFilteredSum <- function(samples) { 
  matrix_list_filtered <- vector(mode = "list", length = 4)
  for (i in 1:length(samples)) { 
    features <- read.delim("features.tsv.gz",header = FALSE, stringsAsFactors = FALSE)
    matrix <- readMM(paste(samples[[i]],"_filtered_matrix.mtx.gz",sep=""))
    barcodes<- read.delim(paste(samples[[i]],"_filtered_barcodes.tsv.gz",sep=""),header = FALSE, stringsAsFactors = FALSE)
    colnames(matrix) <- barcodes$V1
    rownames(matrix) <- features$V2
    matrix <- as.data.frame(as.table(rowSums(matrix)))
    colnames(matrix) <- c("gene_name","counts")
    matrix_list_filtered[[i]] <- matrix
  }
  return(matrix_list_filtered)
}

add_global_label <- function(pwobj, Xlab = NULL, Ylab = NULL, Xgap = 0.03, Ygap = 0.03, ...) {
  ylabgrob <- patchwork::plot_spacer()
  if (!is.null(Ylab)) {
    ylabgrob <- ggplot() +
      geom_text(aes(x = .5, y = .5), label = Ylab, angle = 90, ...) +
      theme_void()
  }
  if (!is.null(Xlab)) {
    xlabgrob <- ggplot() +
      geom_text(aes(x = .5, y = .5), label = Xlab, ...) +
      theme_void()
  }
  if (!is.null(Ylab) & is.null(Xlab)) {
    return((ylabgrob + patchworkGrob(pwobj)) + 
             patchwork::plot_layout(widths = 100 * c(Ygap, 1 - Ygap)))
  }
  if (is.null(Ylab) & !is.null(Xlab)) {
    return((ylabgrob + pwobj) + 
             (xlabgrob) +
             patchwork::plot_layout(heights = 100 * c(1 - Xgap, Xgap),
                                    widths = c(0, 100),
                                    design = "
                                   AB
                                   CC
                                   "
             ))
  }
  if (!is.null(Ylab) & !is.null(Xlab)) {
    return((ylabgrob + pwobj) + 
             (xlabgrob) +
             patchwork::plot_layout(heights = 100 * c(1 - Xgap, Xgap),
                                    widths = 100 * c(Ygap, 1 - Ygap),
                                    design = "
                                   AB
                                   CC
                                   "
             ))
  }
  return(pwobj)
}

renameClusters <- function(seurat.obj,
                           cluster.names,
                           rename.type = c("ident", "annotation.name"),
                           annotation.name = NULL,
                           new.annotation.name = NULL,
                           pull.from = NULL) {
  if (!is.null(pull.from)) {
    current.idents <-
      seurat.obj[[pull.from]] %>%
      as_tibble(rownames = "Cells") %>%
      set_colnames(c("Cells", "Idents"))
  } else {
    current.idents <-
      tibble(Cells = names(Idents(seurat.obj)),
             Idents = Idents(seurat.obj))
  }
  
  annotated.ident <- tibble(Idents = levels(current.idents$Idents),
                            Annotations = cluster.names) %>%
    full_join(current.idents) %>%
    mutate(Annotations = factor(Annotations, levels = unique(Annotations))) %>%
    dplyr::select(-Idents) %>%
    relocate(Cells, .before = Annotations) %>%
    deframe()
  
  if (any(rename.type %in% "ident")) {
    Idents(seurat.obj) <- annotated.ident
  }
  
  if (any(rename.type %in% "annotation.name")) {
    seurat.obj[[annotation.name]] <- annotated.ident
    
  }
  
  if (any(rename.type %in% "new.annotation.name")) {
    seurat.obj[[new.annotation.name]] <- annotated.ident
    
  }
  return(seurat.obj)
  
}

umapCellAnno <- function(seurat.obj,
                         annotation.name = NULL,
                         point.size = 1,
                         label.size = 10,
                         title = "",
                         legend.title.size = 0,
                         legend.text.size = 15,
                         axis.title.x.size = 15,
                         axis.title.y.size = 15,
                         axis.text.y.left.size = 15,
                         axis.text.x.bottom.size = 15,
                         counts.as.title = FALSE,
                         legend = TRUE,
                         cell.legend.size = 10,
                         counts.in.legend = TRUE,
                         use.cols = NULL) {
  if (!is.null(annotation.name)) {
    vars <- annotation.name
  } else if (!class(try(seurat.obj$Seurat_Assignment, silent = TRUE)
  ) == "try-error") {
    vars <- "Seurat_Assignment"
  } else {
    vars <- "seurat_clusters"
  }
  
  umap <-
    as_tibble(Embeddings(seurat.obj, reduction = "umap"), rownames = "Cell") %>%
    mutate(Clusters = FetchData(seurat.obj, vars = vars)[[1]])
  
  cluster.counts <- umap %>%
    group_by_at(4) %>%
    tally() %>%
    arrange(desc(n))
  
  umap$Clusters <-
    factor(umap$Clusters , levels = cluster.counts[[1]])
  umap <- umap[order(umap$Clusters),]
  
  if (is.null(use.cols)) {
    use.cols <-
      hcl(h = seq(15, 375, length = length(unique(umap[[4]])) + 1),
          c = 100,
          l = 65)[seq_along(unique(umap[[4]]))]
  }
  
  if (counts.as.title == TRUE) {
    title = paste(comma(sum(cluster.counts$n)), "Cells")
  }
  
  labels <- cluster.counts[[1]]
  if (isTRUE(counts.in.legend)) {
    labels <-
      as.character(glue("{cluster.counts[[1]]}\n ({cluster.counts[[2]]} Cells)"))
    labels <- paste0(labels, "\n")
  }
  
  l.coord <-
    umap %>% group_by(Clusters) %>% summarize(UMAP1 = median(UMAP_1),
                                              UMAP2 = median(UMAP_2))
  
  p1 <- ggplot(data = umap, mapping = aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = Clusters), size = point.size) +
    scale_color_manual(values = use.cols, labels = labels) +
    theme_bw() +
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank()
    ) +
    theme(
      axis.line = element_line(color = 'black'),
      legend.title = element_text(size = legend.title.size),
      legend.text = element_text(size = legend.text.size),
      axis.title.x = element_text(size = axis.title.x.size),
      axis.title.y = element_text(size = axis.title.y.size),
      axis.text.y.left = element_text(size = axis.text.y.left.size),
      axis.text.x.bottom = element_text(size = axis.text.x.bottom.size)
    ) +
    guides(colour = guide_legend(override.aes = list(size = cell.legend.size))) +
    geom_text_repel(
      data = l.coord,
      mapping = aes(x = UMAP1, y = UMAP2, label = Clusters),
      size = label.size,
      direction = "y"
    )
  return(p1)
}
