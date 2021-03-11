pipeline.createMetacell <- function(seuratObject)
{
  util.info("Creating meta cell")
  patient.cluster <- table( seuratObject$seurat_clusters, seuratObject$orig.ident)
  patient.cluster <- patient.cluster / rowSums(patient.cluster)
  drop.cluster <- names( which( apply( patient.cluster, 1, function(x) sum( sort(x,decreasing=T)[1:2] ) ) > .8 ) )
  # remove patients where xx% of cells belong to one cluster
  # cluster.patient <- table( seuratObject$orig.ident, seuratObject$seurat_clusters )
  # cluster.patient <- cluster.patient[ , which(!colnames(cluster.patient) %in% drop.cluster ) ]
  # cluster.patient <- cluster.patient / rowSums(cluster.patient)
  # barplot(t(cluster.patient))

  labels <- paste0( "c", seuratObject$seurat_clusters, " ", seuratObject$orig.ident )
  names(labels) <- colnames(seuratObject)
  
  # remove meta-cells with less than 10 cells or where 80% of cells belong to one cluster
  labels <- labels[ which(!seuratObject$seurat_clusters%in%drop.cluster) ]
  labels <- labels[ which( labels %in% names( which( table(labels)>=10 ) ) ) ]
  
  sort( table(labels) )
  
  labels.clusterNo <- ceiling( sort( table(labels) ) / 100 )
  #TODO was passiert hier?!
  o <- order( sapply( strsplit(names(labels.clusterNo), "ac" ), function(x) as.numeric(x[2]) )  )
  labels.clusterNo <- labels.clusterNo[o]
  o <- order( sapply( strsplit(names(labels.clusterNo), " " ), function(x) as.numeric(substr(x[1],2,nchar(x[1]))) )  )
  labels.clusterNo <- labels.clusterNo[o]
  
  sum(labels.clusterNo)
  
  metacell.labels <- rep(NA,ncol(seuratObject)) 
  names(metacell.labels) <- colnames(seuratObject)
  metacell.labels[names(labels)] <- labels
  
  seuratObject[["metacellLabelsLvl1"]] <- metacell.labels
  seuratObject[["metacellLabelsLvl2"]] <- metacell.labels
  
  seuratObject[["cellInMetacell"]] <- !is.na(metacell.labels)
  sum(seuratObject$cellInMetacell)
  
  metacell.data <- matrix(NA,nrow(seuratObject),0,dimnames=list(rownames(seuratObject),c()))
  
  #TODO was passiert hier?! achso, x =... ist Testfall für schleife
  x = "c0 ac01"
  pb <-txtProgressBar(min = 0, max = length(labels.clusterNo),style=3)
  
  # create metacell.data, that is the center of kmeans cluster
  for( x in names(labels.clusterNo) )
  {
    mc.cells <- names(which(metacell.labels==x))
    expr <- data.matrix( seuratObject@assays$RNA@data[ , mc.cells ] )
    
    km <- kmeans(t(expr),centers=labels.clusterNo[x])
    
    lab <- paste0( x, " x", seq(max(km$cluster)) )
    
    expr <- t(km$centers)
    colnames(expr) <- lab
    
    metacell.data <- cbind(metacell.data, expr)
    
    seuratObject$metacellLabelsLvl2[mc.cells] <- paste0( seuratObject$metacellLabelsLvl1[mc.cells], " x", km$cluster )
    
    setTxtProgressBar( pb, pb$getVal()+1 )
  }
  pb$kill()
  
  seuratObject[["metacellLabelsLvl1"]] <- as.factor(seuratObject$metacellLabelsLvl1)
  seuratObject[["metacellLabelsLvl2"]] <- as.factor(seuratObject$metacellLabelsLvl2)
  
  return(seuratObject)
}
