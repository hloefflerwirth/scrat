pipeline.summarySheetSeurat <- function(env){
  
  dirname <- paste(env$files.name, "- Results/Summary Sheets - Seurat")
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  
  dirname <- file.path( paste(env$files.name, "- Results"),"Summary Sheets - Seurat","Seurat Clustering.pdf" )
  #TODO
  #util.info("Writing:", file.path(dirname, "*.pdf"))
  #dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  
  pdf(dirname, 21/2.54, 29.7/2.54, useDingbats=FALSE)
  
  
  # Gennamen:
  #env$seuratObject@assays$RNA@var.features
  # numeric metadata
  metadata_numeric =  env$seuratObject@meta.data[, unlist(lapply(env$seuratObject@meta.data, is.numeric))]
  metadata_numeric_names =  colnames(metadata_numeric)
  # factor metadata
  #metadata_factor = env$seuratObject@meta.data[, unlist(lapply(env$seuratObject@meta.data, is.factor))]
  #metadata_factor_names = colnames(metadata_factor)
  metadata_factor_names = c("group.labels", "orig.ident", "seurat_clusters", "metacellLabelsLvl1", "metacellLabelsLvl2")
  metadata_factor = env$seuratObject@meta.data[, metadata_factor_names]
  # get number of levels per factor
  n_levels = mapply(metadata_factor, FUN = function(X){length(levels(X))})

  
  #FeatureScatter(env$seuratObject, feature1 = "nCount_RNA", feature2 = "percent_mt")
  # Plot Metadata with more than one level
  plot_list<- list()
  plot_list[[1]] <- DimPlot(env$seuratObject, reduction = "tsne", group.by = c(metadata_factor_names[(n_levels>1)]), ncol = 3, combine = F)
  plot_list[[2]] <- DimPlot(env$seuratObject, reduction = "pca", group.by = c(metadata_factor_names[(n_levels>1)]), ncol = 3, combine = F)
  plot_list[[3]] <- DimPlot(env$seuratObject, reduction = "umap", group.by = c(metadata_factor_names[(n_levels>1)]), ncol = 3, combine = F)
  
  grid.arrange(grobs = unlist(plot_list, recursive=FALSE), ncol=length(metadata_factor_names[(n_levels>1)]))
  
  plot <- suppressWarnings(FeaturePlot(env$seuratObject, features = c(metadata_numeric_names), reduction = env$preferences$dim.reduction, ncol = 3))
  print(plot)
  
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  seurat.markers <- FindAllMarkers( env$seuratObject, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

  n_genes = 3
  top_genes = seurat.markers %>% group_by(cluster) %>% top_n(n = n_genes, wt = avg_logFC)
  n_cluster = length(levels(top_genes$cluster))

  plot_list<- list()
  #layout(matrix(c(1:(n_genes*n_cluster)), n_cluster, n_genes, byrow =TRUE))
  for(i in as.integer(levels(top_genes$cluster))){
    feature = c(unlist(top_genes["gene"][((i)*n_genes + 1):((i+1)*n_genes),]))
    plot_list[[(i+1)]] <- FeaturePlot(env$seuratObject, features = feature, reduction = env$preferences$dim.reduction, combine = F, ncol = n_genes)
  }
  grid.arrange(grobs = unlist(plot_list, recursive=FALSE), ncol=3)
  #do.call("grid.arrange", c(plot_list, ncol=3))
  dev.off()
}



