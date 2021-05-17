pipeline.createMetacell <- function(env)
{
  util.info("Creating meta cell")
  patient.cluster <- table( env$seuratObject$seurat_clusters, env$seuratObject$orig.ident)
  patient.cluster <- patient.cluster / rowSums(patient.cluster)
  drop.cluster <- names( which( apply( patient.cluster, 1, function(x) sum( sort(x,decreasing=T)[1:2] ) ) > .8 ) )

  labels <- paste0( "c", env$seuratObject$seurat_clusters, " ", env$seuratObject$orig.ident )
  names(labels) <- colnames(env$seuratObject)
  
  # remove meta-cells with less than 10 cells or where 80% of cells belong to one cluster
  labels <- labels[ which(!env$seuratObject$seurat_clusters%in%drop.cluster) ]
  labels <- labels[ which( labels %in% names( which( table(labels)>=10 ) ) ) ]
  
  sort( table(labels) )
  
  labels.clusterNo <- ceiling( sort( table(labels) ) / 100 )
  o <- order( sapply( strsplit(names(labels.clusterNo), "ac" ), function(x) as.numeric(x[2]) )  )
  labels.clusterNo <- labels.clusterNo[o]
  o <- order( sapply( strsplit(names(labels.clusterNo), " " ), function(x) as.numeric(substr(x[1],2,nchar(x[1]))) )  )
  labels.clusterNo <- labels.clusterNo[o]
  
  metacell.labels <- rep(NA,ncol(env$seuratObject)) 
  names(metacell.labels) <- colnames(env$seuratObject)
  metacell.labels[names(labels)] <- labels
  
  env$seuratObject[["metacellLabelsLvl1"]] <- metacell.labels
  env$seuratObject[["metacellLabelsLvl2"]] <- metacell.labels
  
  env$seuratObject[["cellInMetacell"]] <- !is.na(metacell.labels)
  sum(env$seuratObject$cellInMetacell)
  
  metacell.data <- matrix(NA,nrow(env$seuratObject),0,dimnames=list(rownames(env$seuratObject),c()))
  
  pb <-txtProgressBar(min = 0, max = length(labels.clusterNo),style=3)
  
  # create metacell.data, that is the center of kmeans cluster
  for( x in names(labels.clusterNo) )
  {
    mc.cells <- names(which(metacell.labels==x))
    expr <- data.matrix( env$seuratObject@assays$RNA@data[ , mc.cells ] )
    
    km <- kmeans(t(expr),centers=labels.clusterNo[x])
    
    lab <- paste0( x, " x", seq(max(km$cluster)) )
    
    expr <- t(km$centers)
    colnames(expr) <- lab
    
    metacell.data <- cbind(metacell.data, expr)
    
    env$seuratObject$metacellLabelsLvl2[mc.cells] <- paste0( env$seuratObject$metacellLabelsLvl1[mc.cells], " x", km$cluster )
    
    setTxtProgressBar( pb, pb$getVal()+1 )
  }
  pb$kill()
  
  env$metacellData = metacell.data
  env$metacellLabels = names(labels.clusterNo)
  
  env$seuratObject[["metacellLabelsLvl1"]] <- as.factor(env$seuratObject$metacellLabelsLvl1)
  env$seuratObject[["metacellLabelsLvl2"]] <- as.factor(env$seuratObject$metacellLabelsLvl2)
  
  return(env)
}
