pipeline.cellcycleProcessing <- function(env)
{
  util.info("Classification of cell cycle phase")   
  
  if( any(row.names(env$seuratObject) %in% Seurat::cc.genes.updated.2019$g2m.genes) || any(row.names(env$seuratObject) %in% Seurat::cc.genes.updated.2019$s.genes) ){
    marker <- list()
    marker$S <- Seurat::cc.genes.updated.2019$s.genes
    marker$G2M <- Seurat::cc.genes.updated.2019$g2m.genes
  } else
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

  env$seuratObject <- CellCycleScoring(object = env$seuratObject, g2m.features = marker$G2M, s.features = marker$S)
  
  if( any(is.na(env$seuratObject$S.Score)) || any(is.na(env$seuratObject$G2M.Score)))
  {
    util.warn("Cell cycle classificataion failed. Possibly too few features in data set.")
    return(env)
  }
  
  if( env$preferences$activated.modules$reporting )
  {
    filename <- file.path(paste(env$files.name, "- Results"), "Data Overview", "Cell cycle phase.pdf")
    util.info("Writing:", filename)
    pdf(filename, 21/2.54, 21/2.54)

    plot(env$seuratObject$S.Score, env$seuratObject$G2M.Score, xlab="S score", ylab="G2/M score", pch=16, col=env$seuratObject$group.colors, las=1, main="Scores for S and G2/M phases", xlim=c(-1,1), ylim=c(-1,1) )
    lines(x=c(-1,0),y=c(0,0),lty=2)
    lines(x=c(0,0),y=c(-1,0),lty=2)
    lines(x=c(0,1),y=c(0,1),lty=2)
    text(-0.5,0.5,"G2/M",col=rgb(0,0,0,alpha=0.3), cex=1.5)
    text(-0.5,-0.5,"G1",col=rgb(0,0,0,alpha=0.3), cex=1.5)
    text(0.6,0,"S",col=rgb(0,0,0,alpha=0.3), cex=1.5)
    legend("bottomright", as.character(unique(env$group.labels)), cex=0.5, text.col=env$groupwise.group.colors, bg="white")
    
    dev.off()
  }
  
  if (env$preferences$preprocessing$cellcycle.correction)
  {
    util.info("Correction for cell cycle phase") 
    env$seuratObject <- ScaleData(env$seuratObject, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(env$seuratObject))
  }
  return(env)
}

pipeline.prepareIndata <- function(env)
{
  env$indata.sample.mean <- Matrix::colMeans(env$seuratObject)

  if (env$preferences$preprocessing$sample.quantile.normalization)
  {
    #env$indata <- Quantile.Normalization(env$indata)
    env$seuratObject@assays$RNA@data <- Quantile.Normalization(env$seuratObject@assays$RNA@data)
  }

  #colnames(env$seuratObject) <- make.unique(colnames(env$seuratObject))
  names(env$group.labels) <- make.unique(names(env$group.labels))
  names(env$group.colors) <- make.unique(names(env$group.colors))

  env$indata.gene.mean <- Matrix::rowMeans(env$seuratObject)

  if (env$preferences$preprocessing$feature.centralization)
  {
    #env$indata <- env$indata - env$indata.gene.mean
    env$seuratObject@assays$RNA@data = as(env$seuratObject@assays$RNA@data - env$indata.gene.mean, "dgCMatrix")
  }
  return(env)
}


pipeline.generateSOM <- function(env)
{
  env$som.result <- som.linear.init(env$seuratObject@assays$RNA@data,somSize=env$preferences$dim.1stLvlSom)
  
  # Rotate/Flip First lvl SOMs

  if (env$preferences$rotate.SOM.portraits > 0)
  {
    for (i in 1:env$preferences$rotate.SOM.portraits)
    {
      o <- matrix(c(1:(env$preferences$dim.1stLvlSom^2)), env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom, byrow=TRUE)
      o <- o[rev(1:env$preferences$dim.1stLvlSom),]
      env$som.result <- env$som.result[as.vector(o),]
    }
  }

  if (env$preferences$flip.SOM.portraits)
  {
    o <- matrix(c(1:(env$preferences$dim.1stLvlSom^2)), env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom, byrow=TRUE)
    env$som.result <- env$som.result[as.vector(o),]
  }

  env$som.result <- som.training( env$seuratObject@assays$RNA@data, env$som.result, prolongationFactor = env$preferences$training.extension, verbose = TRUE )
    
  env$metadata <- env$som.result$weightMatrix
  colnames(env$metadata) <- colnames(env$seuratObject)

  env$som.result$weightMatrix <- NULL


  ## set up SOM dependent variables
  
  env$gene.info$coordinates <- apply( env$som.result$node.summary[env$som.result$feature.BMU,c("x","y")], 1, paste, collapse=" x " )
  names(env$gene.info$coordinates) <- rownames(env$seuratObject)
  
  return(env)
}
