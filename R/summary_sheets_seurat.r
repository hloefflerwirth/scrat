pipeline.summarySheetSeurat <- function(env)
{
  ## Cell cycle phases

  util.info("Writing preprocessing reports.")

  dir.create( file.path(paste(env$files.name, "- Results"), "Data Overview"), showWarnings=FALSE)

  pt.cex <- 1.2
  if(ncol(env$seuratObject)<10000) pt.cex <- 2
  if(ncol(env$seuratObject)<1000) pt.cex <- 4


  if( !is( try({ length( env$seuratObject$Phase ) }, silent = T), "try-error" ) )
  {
    filename <- file.path(paste(env$files.name, "- Results"), "Data Overview", "Cell cycle phase 1.png")
    png(filename, 1500, 1000)

    plot(env$seuratObject$S.Score, env$seuratObject$G2M.Score, xlab="S score", ylab="G2/M score", pch=16, col=env$seuratObject$group.colors, las=1, main="Scores for S and G2/M phases", xlim=c(-1,1), ylim=c(-1,1), cex=pt.cex )
    lines(x=c(-1,0),y=c(0,0),lty=2)
    lines(x=c(0,0),y=c(-1,0),lty=2)
    lines(x=c(0,1),y=c(0,1),lty=2)
    text(-0.5,0.5,"G2/M",col=rgb(0,0,0,alpha=0.3), cex=1.5)
    text(-0.5,-0.5,"G1",col=rgb(0,0,0,alpha=0.3), cex=1.5)
    text(0.6,0,"S",col=rgb(0,0,0,alpha=0.3), cex=1.5)
    legend("bottomright", as.character(unique(env$seuratObject$group.labels)), cex=0.5, text.col=unique(env$seuratObject$group.colors), bg="white")

    dev.off()


    filename <- file.path(paste(env$files.name, "- Results"), "Data Overview", "Cell cycle phase 2.png")
    png(filename, 1500, 1000)

    par(mfrow=c(1,2))

    par(mar=c(6,5,4,10))
    barplot(table(env$seuratObject$Phase), xlab="cell cycle phase", ylab="cells")

    if( length(unique(env$seuratObject$group.labels))>1 )
    {
      tab <- t(table(env$seuratObject$Phase, env$seuratObject$group.labels))
      tab <- t(apply(tab,1,function(x)x/sum(x)))
      o <- hclust(dist(tab))

      par(mar=c(6,5,4,1))
      image( 1:3, seq(nrow(tab)), z=t(tab[o$order,]),col=env$color.palette.heatmaps(1000),zlim=c(0,1),
             xlab="cell cycle phase", ylab="", axes=F)
      axis(2,seq(nrow(tab)), rownames(tab), tick=F, las=2, cex.axis=.8)
      axis(1,seq(ncol(tab)), colnames(tab), tick=F, las=1, cex.axis=1)
    }

    dev.off()

  }


  filename <- file.path(paste(env$files.name, "- Results"), "Data Overview", "Groups - tSNE.png")
  png(filename, 1500, 1000)

  d <- DimPlot(env$seuratObject, reduction = "tsne", group.by = "group.labels", combine = TRUE, pt.size=pt.cex, cols=env$groupwise.group.colors )
  print(d)

  dev.off()

  filename <- file.path(paste(env$files.name, "- Results"), "Data Overview", "Groups - UMAP.png")
  png(filename, 1500, 1000)

  d <- DimPlot(env$seuratObject, reduction = "umap", group.by = "group.labels", combine = TRUE, pt.size=pt.cex, cols=env$groupwise.group.colors )
  print(d)

  dev.off()

  filename <- file.path(paste(env$files.name, "- Results"), "Data Overview", "Sample IDs - tSNE.png")
  png(filename, 1500, 1000)

  d <- DimPlot(env$seuratObject, reduction = "tsne", group.by = "orig.ident", combine = TRUE, pt.size=pt.cex)
  print(d)

  dev.off()

  filename <- file.path(paste(env$files.name, "- Results"), "Data Overview", "Sample IDs - UMAP.png")
  png(filename, 1500, 1000)

  d <- DimPlot(env$seuratObject, reduction = "umap", group.by = "orig.ident", combine = TRUE, pt.size=pt.cex)
  print(d)

  dev.off()



  filename <- file.path(paste(env$files.name, "- Results"), "Data Overview", "Seurat clusters - tSNE.png")
  png(filename, 1500, 1000)

  d <- DimPlot(env$seuratObject, reduction = "tsne", group.by = "seurat_clusters", combine = TRUE, pt.size=pt.cex)
  print(d)

  dev.off()

  filename <- file.path(paste(env$files.name, "- Results"), "Data Overview", "Seurat clusters - UMAP.png")
  png(filename, 1500, 1000)

  d <- DimPlot(env$seuratObject, reduction = "umap", group.by = "seurat_clusters", combine = TRUE, pt.size=pt.cex)
  print(d)

  dev.off()





  # find markers for every cluster compared to all remaining cells, report only the positive ones
  seurat.markers <- FindAllMarkers( env$seuratObject, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

  if( length(seurat.markers) > 0 )
  {

    n.genes <- 5

    seurat.markers.filter <- by( seurat.markers, seurat.markers$cluster, function(x)
      x <- x[order(x$p_val_adj)[1:n.genes],] )

    seurat.markers.filter <- do.call( c,lapply( seurat.markers.filter, function(x) x$gene ) )
	seurat.markers.filter <- unique( as.vector(na.omit(seurat.markers.filter )) )

    d <- FetchData(env$seuratObject, seurat.markers.filter, slot = "data")
    seurat_clusters <- as.numeric( as.vector(env$seuratObject$seurat_clusters) )


    filename <- file.path(paste(env$files.name, "- Results"), "Data Overview", "Seurat clusters - markers.png")
    png(filename, 1500, 1000)

    layout( matrix(c(0,2,3,1),2), widths=c(1,12), heights=c(1,12) )

    par(mar=c(5,6,1,2))
    image( 1:nrow(d), 1:ncol(d), data.matrix(d)[order(seurat_clusters),],
           axes=F, col=env$color.palette.heatmaps(1000), xlab="cells", ylab="" )

    box()
    x.coord <- sapply( sort(unique(seurat_clusters)), function(x)
      mean( which( sort(seurat_clusters) == x ) )
    )
    axis(1,x.coord,sort(unique(seurat_clusters)), tick=FALSE )

    # y.coord <- seq(3,ncol(d),by=n.genes)
    # axis(2,y.coord,sort(unique(seurat_clusters)),las=2)
    axis(2,1:ncol(d), colnames(d), las=2, tick=FALSE )

    cluster.cols <- color.palette.discrete( max(seurat_clusters)+1 )

    par(mar=c(5,1,1,0))
    image(rbind(1:nrow(d)), col=rep(cluster.cols,each=n.genes), axes=FALSE)

    par(mar=c(0,6,1,2))
    image(cbind(1:ncol(d)), col=cluster.cols[ match( sort(seurat_clusters), sort(unique(seurat_clusters)) ) ], axes=FALSE)

    dev.off()

  }

}


