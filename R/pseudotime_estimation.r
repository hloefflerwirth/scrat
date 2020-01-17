pipeline.pseudotimeEstimation <- function()
{
  if( !preferences$pseudotime.estimation$initiator.sample %in% colnames(metadata) )
  {
    util.warn("Initiator sample not found in data. Using fist cell" )
    preferences$pseudotime.estimation$initiator.sample <<- colnames(metadata)[1]
  }
  
  waypoints <- c( preferences$pseudotime.estimation$initiator.sample, setdiff( colnames(metadata), preferences$pseudotime.estimation$initiator.sample )[seq(1,ncol(metadata)-1,length.out=preferences$pseudotime.estimation$n.waypoints-1)] )
  
  cor.matrix <- cor(metadata)
  diag(cor.matrix) <- 0
  
  adj.matrix.kNNG <- matrix( 0, ncol(metadata), ncol(metadata), dimnames=list(colnames(metadata),colnames(metadata)) )
  
  for( i in 1:ncol(cor.matrix) )
  {
    connect.samples <- which( cor.matrix[,i] >= sort(cor.matrix[,i],decreasing=T)[preferences$pseudotime.estimation$k] )
    for( x in connect.samples )
    {
      adj.matrix.kNNG[connect.samples,i] = cor.matrix[connect.samples,i]
      adj.matrix.kNNG[i,connect.samples] = cor.matrix[i,connect.samples]
    }
  }
  
  if( any( is.infinite( shortest.paths(graph.adjacency(adj.matrix.kNNG,weighted=TRUE),colnames(metadata)[1]) ) ) )
  {
    util.warn("Disabling pseudotime estimation due to isolated subgraphs. Increase k.")
    preferences$pseudotime.estimation <<- NULL
    return()
  }
  
  T.list <- list()
  for( ng.i in seq(preferences$pseudotime.estimation$n.iterations) )
  {
    g.IooK <- graph.empty( ncol(metadata), directed = FALSE)
    V(g.IooK)$name <- colnames(metadata)
    
    for( i in 1:ncol(adj.matrix.kNNG) )
    {
      connect.samples <- sample( names( which( adj.matrix.kNNG[,i] > 0 ) ), preferences$pseudotime.estimation$I )
      g.IooK <- add.edges(g.IooK,as.vector(rbind(colnames(metadata)[i],connect.samples)), attr= list(weight=1-adj.matrix.kNNG[connect.samples,i]) ) 
    }
    
    d.i.j <- shortest.paths(g.IooK,waypoints)
    if( any( is.infinite( d.i.j[1,] ) ) ) next
    
    w.i.j <- t( apply(d.i.j,1,function(x) x^2 / sum(x^2) ) )
    t.i.j <- d.i.j
    for( j in rownames(d.i.j) )
      for( i in colnames(d.i.j) )
      {
        if( d.i.j[1,i] < d.i.j[1,j] ) t.i.j[j,i] <- -t.i.j[j,i]
        t.i.j[j,i] <- t.i.j[j,i] + d.i.j[1,j]
      }
    
    traj.it <- sapply( colnames(d.i.j),function(i)  sum(t.i.j[,i]*w.i.j[,i]/preferences$pseudotime.estimation$n.waypoints) )
    it <- 0
    while(TRUE) 
    {
      t.i.j <- d.i.j
      for( j in rownames(d.i.j) )
        for( i in colnames(d.i.j) )
        {
          if( traj.it[i] < traj.it[j] ) t.i.j[j,i] <- -t.i.j[j,i]
          t.i.j[j,i] <- t.i.j[j,i] + d.i.j[1,j]
        }
      
      traj.last <- traj.it
      traj.it <- sapply( colnames(d.i.j),function(i)  sum(t.i.j[,i]*w.i.j[,i]/preferences$pseudotime.estimation$n.waypoints) )
      
      it <- it+1
      if( cor( traj.last, traj.it ) == 1 || 
          cor( traj.last, traj.it ) > 0.95 && it > 100 ||
          it > 200 ) break
      print(it)
    }
    
    T.list[[ng.i]] <- traj.it
  }
  
  n.skipped.elements <- sum( sapply(T.list,is.null) )
  if( n.skipped.elements > 0 )
  {
    util.warn("Skipped",n.skipped.elements,"iteration(s) due to isolated subgraphs.")
  }
  if( length(T.list) == 0 )
  {
    util.warn("Disabling pseudotime estimation. Increase I.")
    preferences$pseudotime.estimation <<- NULL
    return()
  }
  
  pseudotime.trajectory <<- colMeans( do.call(rbind,T.list) )
  pseudotime.trajectory <<- (pseudotime.trajectory-min(pseudotime.trajectory))/(max(pseudotime.trajectory)-min(pseudotime.trajectory))
}


