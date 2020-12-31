pipeline.pseudotimeEstimation <- function(env)
{
  if( !env$preferences$pseudotime.estimation$initiator.sample %in% colnames(env$metadata) )
  {
    util.warn("Initiator sample not found in data. Using fist cell" )
    env$preferences$pseudotime.estimation$initiator.sample <- colnames(env$metadata)[1]
  }
  
  waypoints <- c( env$preferences$pseudotime.estimation$initiator.sample, setdiff( colnames(env$metadata), env$preferences$pseudotime.estimation$initiator.sample )[seq(1,ncol(env$metadata)-1,length.out=env$preferences$pseudotime.estimation$n.waypoints-1)] )
  
  cor.matrix <- cor(env$metadata)
  diag(cor.matrix) <- 0
  
  adj.matrix.kNNG <- matrix( 0, ncol(env$metadata), ncol(env$metadata), dimnames=list(colnames(env$metadata),colnames(env$metadata)) )
  
  for( i in 1:ncol(cor.matrix) )
  {
    connect.samples <- which( cor.matrix[,i] >= sort(cor.matrix[,i],decreasing=T)[env$preferences$pseudotime.estimation$k] )
    for( x in connect.samples )
    {
      adj.matrix.kNNG[connect.samples,i] = cor.matrix[connect.samples,i]
      adj.matrix.kNNG[i,connect.samples] = cor.matrix[i,connect.samples]
    }
  }
  
  if( any( is.infinite( shortest.paths(graph.adjacency(adj.matrix.kNNG,weighted=TRUE),colnames(env$metadata)[1]) ) ) )
  {
    util.warn("Disabling pseudotime estimation due to isolated subgraphs. Increase k.")
    env$preferences$pseudotime.estimation <- NULL
    return(env)
  }
  
  T.list <- list()
  for( ng.i in seq(env$preferences$pseudotime.estimation$n.iterations) )
  {
    g.IooK <- graph.empty( ncol(env$metadata), directed = FALSE)
    V(g.IooK)$name <- colnames(env$metadata)
    
    for( i in 1:ncol(adj.matrix.kNNG) )
    {
      connect.samples <- sample( names( which( adj.matrix.kNNG[,i] > 0 ) ), env$preferences$pseudotime.estimation$I )
      g.IooK <- add.edges(g.IooK,as.vector(rbind(colnames(env$metadata)[i],connect.samples)), attr= list(weight=1-adj.matrix.kNNG[connect.samples,i]) ) 
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
    
    traj.it <- sapply( colnames(d.i.j),function(i)  sum(t.i.j[,i]*w.i.j[,i]/env$preferences$pseudotime.estimation$n.waypoints) )
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
      traj.it <- sapply( colnames(d.i.j),function(i)  sum(t.i.j[,i]*w.i.j[,i]/env$preferences$pseudotime.estimation$n.waypoints) )
      
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
    env$preferences$pseudotime.estimation <- NULL
    return(env)
  }
  
  env$pseudotime.trajectory <- colMeans( do.call(rbind,T.list) )
  env$pseudotime.trajectory <- (pseudotime.trajectory-min(pseudotime.trajectory))/(max(pseudotime.trajectory)-min(pseudotime.trajectory))
  return(env)
}


