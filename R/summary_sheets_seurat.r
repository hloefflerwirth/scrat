pipeline.summarySheetSeurat <- function(env){
  
  filename <- file.path(paste(env$files.name, "- Results"), "Data Overview", "Clustering.pdf")
  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 29.7/2.54, useDingbats=FALSE)
  
  # Gennamen:
  #env$seuratObject@assays$RNA@var.features
  # numeric metadata
  metadata_numeric =  env$seuratObject@meta.data[, unlist(lapply(env$seuratObject@meta.data, is.numeric))]
  metadata_numeric_names =  colnames(metadata_numeric)
  # factor metadata
  #metadata_factor = env$seuratObject@meta.data[, unlist(lapply(env$seuratObject@meta.data, is.factor))]
  #metadata_factor_names = colnames(metadata_factor)
  metadata_factor_names = c("group.labels", "orig.ident", "seurat_clusters")
  metadata_factor = env$seuratObject@meta.data[, metadata_factor_names]
  # get number of levels per factor
  n_levels = mapply(metadata_factor, FUN = function(X){length(levels(X))})

  # Plot Metadata with more than one level
  plot_list<- list()
  plot_list[[1]] <- DimPlot(env$seuratObject, reduction = "tsne", group.by = c(metadata_factor_names[(n_levels>1)]), combine = F)
  plot_list[[2]] <- DimPlot(env$seuratObject, reduction = "umap", group.by = c(metadata_factor_names[(n_levels>1)]),  combine = F)
  
  grid.arrange(grobs = unlist(plot_list, recursive=FALSE), ncol=2)
  
  if(env$preferences$preprocessing$create.meta.cell){
    metadata_factor_names = c("metacellLabelsLvl1", "metacellLabelsLvl2", "cellInMetacell")
    metadata_factor = env$seuratObject@meta.data[, metadata_factor_names]
    # get number of levels per factor
    n_levels = mapply(metadata_factor, FUN = function(X){length(levels(X))})
    
    # Plot Metadata with more than one level
    plot_list<- list()
    plot_list[[1]] <- DimPlot(env$seuratObject, reduction = "tsne", group.by = c(metadata_factor_names[(n_levels>1)]), combine = F)

    grid.arrange(grobs = unlist(plot_list, recursive=FALSE), ncol=1)
  }
  
  # plot of numerical meta data
  plot <- suppressWarnings(FeaturePlot(env$seuratObject, features = c(metadata_numeric_names), reduction = env$preferences$dim.reduction, ncol = 1))
  print(plot)
  
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  seurat.markers <- FindAllMarkers( env$seuratObject, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

  n_genes = 3
  top_genes = seurat.markers %>% group_by(cluster) %>% top_n(n_genes)
  n_cluster = length(levels(top_genes$cluster))

  plot_list<- list()
  #layout(matrix(c(1:(n_genes*n_cluster)), n_cluster, n_genes, byrow =TRUE))
  for(i in 0:(length(unique(top_genes$cluster))-1)){
    feature = c(unlist(top_genes["gene"][((i)*n_genes + 1):((i+1)*n_genes),]))
    plot_list[[(i+1)]] <- FeaturePlot(env$seuratObject, features = feature, reduction = env$preferences$dim.reduction, combine = F, ncol = n_genes)
  }
  #label <- textGrob(c("", "line1", "line2"))
  #plots <- arrangeGrob(unlist(plot_list, recursive=FALSE), ncol=3)
  #grid.arrange(grobs = unlist(list(label, plots)),top = "Top three positive marker genes per cluster")
  grid.arrange(grobs = unlist(plot_list, recursive=FALSE), ncol=3,top = "Top three positive marker genes per cluster", left = "Seurat Clusters")
  dev.off()
}


