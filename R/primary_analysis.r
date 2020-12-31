pipeline.countProcessing <- function(env)
{
  if (env$preferences$preprocessing$count.processing)
  {
    util.info("Read count preprocessing")    
    
  }
}


pipeline.cellcycleProcessing <- function(env)
{
  util.info("Classification of cell cycle phase")   
  
  if( length(grep("ENSMUSG",env$gene.info$ids)>0) )
  {
    marker.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))    
  } else
  if( length(grep("ENSG",env$gene.info$ids)>0) )
  {
    marker.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))  
  } else
  {
    util.warn("Cell cycle markers only available for human and mouse organisms")
    return(env)
  }

  cc.scores <- cyclone( env$indata[ which( env$gene.info$ids != "" ), ], marker.pairs, gene.names=env$gene.info$ids[ which( env$gene.info$ids != "" ) ] )
  
  if( any(is.na(cc.scores$scores)) )
  {
    util.warn("Cell cycle classificataion failed. Possibly too few features in data set.")
    return(env)
  }
  
  if( env$preferences$activated.modules$reporting )
  {
    filename <- file.path(paste(env$files.name, "- Results"), "Data Overview", "Cell cycle phase.pdf")
    util.info("Writing:", filename)
    pdf(filename, 21/2.54, 21/2.54)
  
    plot(cc.scores$score$G1, cc.scores$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16, col=env$group.colors, las=1, main="Scores for G1 and G2/M phases", xlim=c(0,1), ylim=c(0,1) )
      abline(h=0.5,lty=2)
      abline(v=0.5,lty=2)
      text(0.75,0.25,"G1",col=rgb(0,0,0,alpha=0.3), cex=1.5)
      text(0.25,0.75,"G2/M",col=rgb(0,0,0,alpha=0.3), cex=1.5)
      text(0.25,0.25,"S",col=rgb(0,0,0,alpha=0.3), cex=1.5)
    
    dev.off()
  }
  
  if (env$preferences$preprocessing$cellcycle.correction)
  {
    util.info("Correction for cell cycle phase") 
    
    env$indata <- removeBatchEffect(env$indata, covariates=cc.scores$score[,c("G1", "G2M")])
  }
  return(env)
}

pipeline.prepareIndata <- function(env)
{
  env$indata.sample.mean <- colMeans(env$indata)

  if (env$preferences$preprocessing$sample.quantile.normalization)
  {
    env$indata <- Quantile.Normalization(env$indata)
  }

  colnames(env$indata) <- make.unique(colnames(env$indata))
  names(env$group.labels) <- make.unique(names(env$group.labels))
  names(env$group.colors) <- make.unique(names(env$group.colors))


  env$indata.gene.mean <- rowMeans(env$indata)

  if (env$preferences$preprocessing$feature.centralization)
  {
    env$indata <- env$indata - env$indata.gene.mean
  }
  return(env)
}


pipeline.generateSOM <- function(env)
{
  env$som.result <- som.linear.init(env$indata,somSize=env$preferences$dim.1stLvlSom)
  
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

  env$som.result <- som.training( env$indata, env$som.result, prolongationFactor = env$preferences$training.extension, verbose = TRUE )
    
  env$metadata <- env$som.result$weightMatrix
  colnames(env$metadata) <- colnames(env$indata)

  env$som.result$weightMatrix <- NULL


  ## set up SOM dependent variables
  
  env$gene.info$coordinates <- apply( env$som.result$node.summary[env$som.result$feature.BMU,c("x","y")], 1, paste, collapse=" x " )
  names(env$gene.info$coordinates) <- rownames(env$indata)
  
  return(env)
}
