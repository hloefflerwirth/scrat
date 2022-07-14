pipeline.seuratPreprocessing <- function(env)
{
  env$seuratObject[["percent.mt"]] <- PercentageFeatureSet(object = env$seuratObject, assay = "RNA", pattern = "^MT-")
    
  # Filtering
  # 1) Low-quality cells where more than 15% of the read counts come from mitochondrial genome
  # 2) Cells with more than 7000 detected genes
  
  dropped.cols = which(env$seuratObject$percent.mt >= 15 | env$seuratObject$nFeature_RNA >= 7000)

  env$seuratObject <- subset(env$seuratObject, subset = nFeature_RNA < 7000 & percent.mt < 15)
    
  if (length(dropped.cols) > 0)
  {
    env$group.labels <- env$group.labels[-dropped.cols]
    env$group.colors <- env$group.colors[-dropped.cols]
    util.warn("Filtered",length(dropped.cols),"cells from data set in preprocessing.")
  }
    
  # Normalization & reduction
  env$seuratObject <- NormalizeData(env$seuratObject, assay = 'RNA', normalization.method = "LogNormalize")
  env$seuratObject <- FindVariableFeatures(env$seuratObject, selection.method = "vst", assay = 'RNA')
  env$seuratObject <- ScaleData(env$seuratObject, model.use = 'linear', vars.to.regress = 'orig.ident')
  
  env$seuratObject <- RunPCA(env$seuratObject, npcs = min(ncol(env$seuratObject)-1, 100))
  env$seuratObject <- RunTSNE(env$seuratObject, reduction = "pca", assay = 'RNA', dims = 1:100, perplexity = 5)
  env$seuratObject <- RunUMAP(env$seuratObject, reduction = "pca", assay = 'RNA', dims = 1:100)
  
  # primary cell clustering
  env$seuratObject <- FindNeighbors(env$seuratObject, assay = 'RNA',  dims = 1:100)
  env$seuratObject <- FindClusters(env$seuratObject, resolution = 1)
  
  return(env)
}


pipeline.cellcycleProcessing <- function(env)
{
  if( any(row.names(env$seuratObject) %in% Seurat::cc.genes.updated.2019$g2m.genes) && any(row.names(env$seuratObject) %in% Seurat::cc.genes.updated.2019$s.genes) ){
    marker <- list()
    marker$S <- intersect( Seurat::cc.genes.updated.2019$s.genes, rownames(env$seuratObject) )
    marker$G2M <- intersect( Seurat::cc.genes.updated.2019$g2m.genes, rownames(env$seuratObject) )
  } 
  else
  {
    if (!biomart.available(env))
    {
      util.warn("Requested biomaRt host seems to be down.")
      util.warn("Disabling classification of cell cycle phase.")
      return(env)
    }
    
    if( length(grep("ENSMUSG",row.names(env$seuratObject))>0) )
    {
      # convert human marker names to mouse gene ids
      convertHumanGeneList <- function(x){
        human = useMart(biomart=env$preferences$database.biomart, dataset = "hsapiens_gene_ensembl", host=env$preferences$database.host)
        mouse = useMart(biomart=env$preferences$database.biomart, dataset = "mmusculus_gene_ensembl", host=env$preferences$database.host)
        genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("ensembl_gene_id"), martL = mouse, uniqueRows=T)
        humanx <- unique(genesV2[, 2])
        return(humanx)
      }
      marker <- list()
      marker$S <- convertHumanGeneList(Seurat::cc.genes.updated.2019$s.genes)
      marker$G2M <- convertHumanGeneList(Seurat::cc.genes.updated.2019$g2m.genes)
    } else
      if( length(grep("ENSG",row.names(env$seuratObject))>0) )
      {
        # convert marker names to gene ids
        convertGeneList <- function(x){
          human = useMart(biomart=env$preferences$database.biomart, dataset = "hsapiens_gene_ensembl", host=env$preferences$database.host)
          genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=T)
          gene_ids <- unique(genesV2[, 2])
          return(gene_ids)
        }
        marker <- list()
        marker$S <- convertGeneList(Seurat::cc.genes.updated.2019$s.genes)
        marker$G2M <- convertGeneList(Seurat::cc.genes.updated.2019$g2m.genes)
      } else
      {
        util.warn("Cell cycle markers only available for human and mouse organisms")
        return(env)
      }
  }
  
  env$seuratObject <- CellCycleScoring(object = env$seuratObject, g2m.features = marker$G2M, s.features = marker$S)
  
  if( any(is.na(env$seuratObject$S.Score)) || any(is.na(env$seuratObject$G2M.Score)))
  {
    util.warn("Cell cycle classificataion failed. Possibly too few features in data set.")
    return(env)
  }
  
  if (env$preferences$preprocessing$cellcycle.correction)
  {
    util.info("Correction for cell cycle phase") 
    env$seuratObject <- ScaleData(env$seuratObject, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(env$seuratObject))
  }
  return(env)
}
