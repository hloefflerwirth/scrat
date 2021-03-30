pipeline.seuratPreprocessing <- function(env)
{
  util.info("Creating and preprocessing Seurat Object")
  colnames(env$indata) <- make.unique(colnames(env$indata))
  env$seuratObject  <- CreateSeuratObject(env$indata, project = env$preferences$dataset.name)
  remove(indata, envir = env)
  # preprocess count matrix
  if(env$preferences$indata.counts == TRUE){  
    # Quality-Control der counts 
    env$seuratObject[["percent.mt"]] <- PercentageFeatureSet(env$seuratObject, pattern = "^MT-")
    
    # 1) Low-quality cells where more than 15% of the read counts derived from the mitochondrial genome were also discarded.
    # 2) Cells with more than 7000 detected genes were discarded to remove likely doublet or multiplet captures.
    dropped.cols = which(env$seuratObject$percent.mt >= 15 | env$seuratObject$nFeature_RNA >= 7000)
    env$seuratObject <- subset(env$seuratObject, subset = nFeature_RNA < 7000 & percent.mt < 15)
    
    if (length(dropped.cols) > 0)
    {
      env$group.labels <- env$group.labels[-dropped.cols]
      env$group.colors <- env$group.colors[-dropped.cols]
      util.warn("Removed",length(dropped.cols),"columns from data set in preprocessing.")
    }
    
    # Normalisierung der Counts
    if (env$preferences$preprocessing$seurat.normalize)
    {
      env$seuratObject <- NormalizeData(env$seuratObject, assay = 'RNA')
    }
  }
  env$seuratObject[["group.labels"]] <- as.factor(env$group.labels)
  env$seuratObject[["group.colors"]] <- env$group.colors
  #reduktion
  env$seuratObject <- FindVariableFeatures(env$seuratObject, assay = 'RNA')
  env$seuratObject <- ScaleData(env$seuratObject)
  env$seuratObject <- RunPCA(env$seuratObject, npcs = min(dim(env$seuratObject)[2]-1, 100))
  env$seuratObject <- RunTSNE(env$seuratObject, reduction = "pca", assay = 'RNA', dims = 1:50, perplexity = 5)
  env$seuratObject <- RunUMAP(env$seuratObject, assay = 'RNA', dims = 1:50)
  
  # cluster von ähnlichen zellen
  env$seuratObject <- FindNeighbors(env$seuratObject, assay = 'RNA',  dims = 1:50)
  env$seuratObject <- FindClusters(env$seuratObject, resolution = 1)
  
  if (env$preferences$preprocessing$create.meta.cell) {
    env$seuratObject <- pipeline.createMetacell(env$seuratObject)
  }
  return(env)
}