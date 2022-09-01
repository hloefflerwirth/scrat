pipeline.createMetacells <- function(env)
{
  util.info("Creating meta cell")
  patient.cluster <- table( env$seuratObject$seurat_clusters, env$seuratObject$orig.ident)
  patient.cluster <- patient.cluster / rowSums(patient.cluster)

  labels <- paste0( "c", env$seuratObject$seurat_clusters, " ", env$seuratObject$orig.ident )
  names(labels) <- colnames(env$seuratObject)

  # remove cluster where >80% of cells come from one or two samples/patient
      # drop.cluster <- names( which( apply( patient.cluster, 1, function(x) sum( sort(x,decreasing=T)[1:2] ) ) > .8 ) )
      # labels <- labels[ which(!env$seuratObject$seurat_clusters%in%drop.cluster) ]

  # remove meta-cells with less than 10 cells

  if( any(table(labels)>=10) )
    labels <- labels[ which( labels %in% names( which( table(labels)>=10 ) ) ) ]


  if( length(unique(labels))==1 )
  {
    util.warn("Meta-cell computation would lead to only one meta-cell. Skipped.")
    return(env)
  }


  labels.clusterNo <- ceiling( sort( table(labels) ) / 100 )
  o <- order( sapply( strsplit(names(labels.clusterNo), " " ), function(x) x[2] )  )
  labels.clusterNo <- labels.clusterNo[o]
  o <- order( sapply( strsplit(names(labels.clusterNo), " " ), function(x) as.numeric(substr(x[1],2,nchar(x[1]))) )  )
  labels.clusterNo <- labels.clusterNo[o]

  metacell.labels.lvl1 <- rep(NA,ncol(env$seuratObject))
  names(metacell.labels.lvl1) <- colnames(env$seuratObject)
  metacell.labels.lvl1[names(labels)] <- labels

  env$seuratObject[["metacell.labels"]] <- metacell.labels.lvl1
  env$seuratObject[["cellInMetacell"]] <- !is.na(metacell.labels.lvl1)

  metacell.data <- matrix(NA,nrow(env$seuratObject),0,dimnames=list(rownames(env$seuratObject),c()))

  pb <-txtProgressBar(min = 0, max = length(labels.clusterNo),style=3)

  # create metacell.data, that is the center of kmeans cluster
  for( x in names(labels.clusterNo) )
  {
    mc.cells <- names(which(metacell.labels.lvl1==x))
    expr <- data.matrix( GetAssayData( env$seuratObject[ , mc.cells ], slot = "data") )

    km <- kmeans(t(expr),centers=labels.clusterNo[x])

    expr <- t(km$centers)
    colnames(expr) <- paste0( x, " x", seq(max(km$cluster)) )

    metacell.data <- cbind(metacell.data, expr)

    env$seuratObject$metacell.labels[mc.cells] <- paste0( metacell.labels.lvl1[mc.cells], " p", km$cluster )

    setTxtProgressBar( pb, pb$getVal()+1 )
  }
  pb$kill()

  env$metacell.data <- metacell.data
  env$group.labels <- rep("auto",ncol(env$metacell.data))
  names(env$group.labels) <- colnames(env$metacell.data)
  env$group.colors <- rep("#000000", ncol(env$metacell.data))
  names(env$group.colors) <- colnames(env$metacell.data)

  util.info(paste(ncol(env$metacell.data),"meta-cells were computed.")

  return(env)
}
