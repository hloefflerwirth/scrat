pipeline.countProcessing <- function()
{
  if (preferences$preprocessing$count.processing)
  {
    util.info("Read count preprocessing")    
    
  #   sce <- SingleCellExperiment( assays=list(counts=indata) )
  #   sce <- calculateQCMetrics(sce)
  # 
  #   libsize.drop <- isOutlier(sce$total_counts, n=3, type="lower", log=TRUE)
  #   feature.drop <- isOutlier(sce$total_features, n=3, type="lower")
  #   expression.drop <- rowMeans(indata) < .1
  # 
  #   util.info("Filtering:",sum(libsize.drop|feature.drop),"cells and",sum(expression.drop),"genes removed.")
  # 
  #   if( env$preferences$activated.modules$reporting )
  #   {
  #     dirname <- paste(files.name, "- Results/Data Overview")
  #     dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  # 
  #     filename <- file.path(paste(files.name, "- Results"), "Data Overview", "Read count preprocessing.pdf")
  #     util.info("Writing:", filename)
  #     pdf(filename, 29.7/2.54, 21/2.54)
  # 
  #     layout(matrix(c(1,3,2,1,5,4),3),heights=c(1,1,15))
  #     par(mar=c(0,0,0,0)); frame()
  #     text(0.5,0.5,"Histograms of library sizes and number of expressed genes",cex=2.5)
  # 
  #     par(mar=c(5,5,0,1))
  #     hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="", las=1, cex.axis=1.4,
  #          breaks=20, col="grey80", ylab="Number of cells", cex.lab=1.4)
  #       abline( v= min( sce$total_counts[!libsize.drop] )/1e6, lty=2, col="red3" )
  #     par(mar=c(0,0,0,1)); frame()
  #       legend("bottomright", "libsize filter threshold\n(3*MAD criterion)", lty=2, col="red3")
  # 
  #     par(mar=c(5,5,0,1))
  #     hist(sce$total_features, xlab="Number of expressed genes", main="", las=1, cex.axis=1.4,
  #        breaks=20, col="grey80", ylab="Number of cells", cex.lab=1.4)
  #       abline( v= min( sce$total_features[!feature.drop] ), lty=2, col="red3" )
  #     par(mar=c(0,0,0,1)); frame()
  #       legend("bottomright", "gene number filter threshold\n(3*MAD criterion)", lty=2, col="red3")
  #   }
  # 
  #   sce <- sce[, which(!libsize.drop&!feature.drop) ]
  #   sce <- sce[ which(!expression.drop), ]
  # 
  #   group.labels <<- group.labels[which(!libsize.drop&!feature.drop)]
  #   group.colors <<- group.colors[which(!libsize.drop&!feature.drop)]
  # 
  #   util.info("Scale Factor normalization")
  # 
  #   suppressWarnings({ sce <- computeSumFactors(sce, sizes=seq(1,ncol(sce),length.out=4)) })
  # 
  #   if( env$preferences$activated.modules$reporting )
  #   {
  #     layout(matrix(c(1,0,1,2,1,0),2),widths=c(1,5,1),heights=c(1,16) )
  # 
  #     par(mar=c(0,0,0,0)); frame()
  #     text(0.5,0.5,"Normalization using scale factors",cex=2.5)
  # 
  #     par(mar=c(5,6,4,5))
  #     plot(sizeFactors(sce), sce$total_counts/1e6, log="xy", ylab="Library size (millions)", xlab="Scale (size) factor", pch=16, col=group.colors, cex.lab=1.4, cex.axis=1.4,las=1)
  # 
  #     dev.off()
  #   }
  # 
  #   sce <- normalize(sce)
  #   dn <- dimnames(sce@assays[["counts"]])
  #   indata <<- sce@assays[["logcounts"]]
  #   dimnames(indata) <<- dn
  }
}


pipeline.cellcycleProcessing <- function()
{
  util.info("Classification of cell cycle phase")   
  
  if( length(grep("ENSMUSG",gene.info$ids)>0) )
  {
    marker.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))    
  } else
  if( length(grep("ENSG",gene.info$ids)>0) )
  {
    marker.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))  
  } else
  {
    util.warn("Cell cycle markers only available for human and mouse organisms")
    return()
  }

  cc.scores <- cyclone( indata[ which( gene.info$ids != "" ), ], marker.pairs, gene.names=gene.info$ids[ which( gene.info$ids != "" ) ] )
  
  if( any(is.na(cc.scores$scores)) )
  {
    util.warn("Cell cycle classificataion failed. Possibly too few features in data set.")
    return()
  }
  
  if( env$preferences$activated.modules$reporting )
  {
    filename <- file.path(paste(files.name, "- Results"), "Data Overview", "Cell cycle phase.pdf")
    util.info("Writing:", filename)
    pdf(filename, 21/2.54, 21/2.54)
  
    plot(cc.scores$score$G1, cc.scores$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16, col=group.colors, las=1, main="Scores for G1 and G2/M phases", xlim=c(0,1), ylim=c(0,1) )
      abline(h=0.5,lty=2)
      abline(v=0.5,lty=2)
      text(0.75,0.25,"G1",col=rgb(0,0,0,alpha=0.3), cex=1.5)
      text(0.25,0.75,"G2/M",col=rgb(0,0,0,alpha=0.3), cex=1.5)
      text(0.25,0.25,"S",col=rgb(0,0,0,alpha=0.3), cex=1.5)
    
    dev.off()
  }
  
  if (preferences$preprocessing$cellcycle.correction)
  {
    util.info("Correction for cell cycle phase") 
    
    indata <<- removeBatchEffect(indata, covariates=cc.scores$score[,c("G1", "G2M")])
  }
  
}

pipeline.prepareIndata <- function()
{
  indata.sample.mean <<- colMeans(indata)

  if (preferences$preprocessing$sample.quantile.normalization)
  {
    indata <<- Quantile.Normalization(indata)
  }

  colnames(indata) <<- make.unique(colnames(indata))
  names(group.labels) <<- make.unique(names(group.labels))
  names(group.colors) <<- make.unique(names(group.colors))


  indata.gene.mean <<- rowMeans(indata)

  if (preferences$preprocessing$feature.centralization)
  {
    indata <<- indata - indata.gene.mean
  }
}


pipeline.generateSOM <- function()
{
  som.result <<- som.linear.init(indata,somSize=preferences$dim.1stLvlSom)
  
  # Rotate/Flip First lvl SOMs

  if (preferences$rotate.SOM.portraits > 0)
  {
    for (i in 1:preferences$rotate.SOM.portraits)
    {
      o <- matrix(c(1:(preferences$dim.1stLvlSom^2)), preferences$dim.1stLvlSom, preferences$dim.1stLvlSom, byrow=TRUE)
      o <- o[rev(1:preferences$dim.1stLvlSom),]
      som.result <<- som.result[as.vector(o),]
    }
  }

  if (preferences$flip.SOM.portraits)
  {
    o <- matrix(c(1:(preferences$dim.1stLvlSom^2)), preferences$dim.1stLvlSom, preferences$dim.1stLvlSom, byrow=TRUE)
    som.result <<- som.result[as.vector(o),]
  }

  som.result <<- som.training( indata, som.result, prolongationFactor = preferences$training.extension, verbose = TRUE )
    
  metadata <<- som.result$weightMatrix
  colnames(metadata) <<- colnames(indata)

  som.result$weightMatrix <<- NULL


  ## set up SOM dependent variables
  
  gene.info$coordinates <<- apply( som.result$node.summary[som.result$feature.BMU,c("x","y")], 1, paste, collapse=" x " )
  names(gene.info$coordinates) <<- rownames(indata)
  
}
