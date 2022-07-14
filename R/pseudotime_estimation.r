pipeline.pseudotimeEstimation <- function(env)
{
  if( length(env$preferences$pseudotime.estimation)==1 && env$preferences$pseudotime.estimation == TRUE ) # can be TRUE -> default parameters
  {
    if(ncol(env$metadata)>=10)
    {
      util.info("Using standard parameters for pseudotime estimation.")
      env$preferences$pseudotime.estimation <- list( 
        n.waypoints = max( 5, min( ceiling(ncol(env$metadata)*0.05), 50 ) ),
        n.iterations = 20,
        k = max( 5, min( ceiling(ncol(env$metadata)*0.05), 50 ) ),
        I = max( 2, min( ceiling(ncol(env$metadata)*0.01), 10 ) ),
        initiator.sample = colnames(env$metadata)[1])
    } else
    {
      util.warn("Disabling pseudotime estimation. Too few samples.")
      env$preferences$pseudotime.estimation <- NULL
    }
  } 
  else
  if( is.list(env$preferences$pseudotime.estimation) ) # use user defined parameters
  {
    if(ncol(env$metadata)>=10)
    {
      default <- list( n.waypoints = max( 5, min( ceiling(ncol(env$metadata)*0.05), 50 ) ),
                       n.iterations = 20,
                       k = max( 5, min( ceiling(ncol(env$metadata)*0.05), 50 ) ),
                       I = max( 2, min( ceiling(ncol(env$metadata)*0.01), 10 ) ),
                       initiator.sample = colnames(env$metadata)[1])
      
      env$preferences$pseudotime.estimation <- modifyList( default, env$preferences$pseudotime.estimation )
      
      if (!is.numeric(env$preferences$pseudotime.estimation$n.waypoints) ||
          env$preferences$pseudotime.estimation$n.waypoints < 3 ||
          env$preferences$pseudotime.estimation$n.waypoints > ncol(env$metadata))
      {
        util.warn("Invalid value of \"pseudotime.estimation$n.waypoints\". Using standard value")
        env$preferences$pseudotime.estimation$n.waypoints <- max( 5, min( ceiling(ncol(env$metadata)*0.05), 50 ) )
      }
      if (!is.numeric(env$preferences$pseudotime.estimation$n.iterations) ||
          env$preferences$pseudotime.estimation$n.iterations < 2 )
      {
        util.warn("Invalid value of \"pseudotime.estimation$n.iterations\". Using 20")
        env$preferences$pseudotime.estimation$n.iterations <- 20
      }
      if (!is.numeric(env$preferences$pseudotime.estimation$k) ||
          env$preferences$pseudotime.estimation$k < 3 ||
          env$preferences$pseudotime.estimation$k > ncol(env$metadata))
      {
        util.warn("Invalid value of \"pseudotime.estimation$k\". Using standard value")
        env$preferences$pseudotime.estimation$k <- max( 5, min( ceiling(ncol(env$metadata)*0.05), 50 ) )
      }
      if (!is.numeric(env$preferences$pseudotime.estimation$I) ||
          env$preferences$pseudotime.estimation$I < 2 ||
          env$preferences$pseudotime.estimation$I >= env$preferences$pseudotime.estimation$k)
      {
        util.warn("Invalid value of \"pseudotime.estimation$I\". Using standard value")
        env$preferences$pseudotime.estimation$I <- max( 2, min( ceiling(ncol(env$metadata)*0.01), 10 ) )
      }
      if( is.numeric(env$preferences$pseudotime.estimation$initiator.sample) && length(env$preferences$pseudotime.estimation$initiator.sample) == 1 && env$preferences$pseudotime.estimation$initiator.sample %in% c(1:ncol(env$metadata)) )
      {
        env$preferences$pseudotime.estimation$initiator.sample <- colnames(env$metadata)[env$preferences$pseudotime.estimation$initiator.sample]
      } 
      if( !is.character(env$preferences$pseudotime.estimation$initiator.sample) || length(env$preferences$pseudotime.estimation$initiator.sample) != 1 || !env$preferences$pseudotime.estimation$initiator.sample %in% colnames(env$metadata) )
      {
        util.warn("Invalid value of \"pseudotime.estimation$initiator.sample\". Using fist sample." )
        env$preferences$pseudotime.estimation$initiator.sample <- colnames(env$metadata)[1]
      }        
    } else
    {
      util.warn("Disabling pseudotime estimation. Too few samples.")
      env$preferences$pseudotime.estimation <- NULL
    }    
  } else
  {
    env$preferences$pseudotime.estimation <- NULL
  } 

  
  
  if( !is.null(env$metacell.data) )
  {
    util.info("Using fist meta-cell as root." )
    env$preferences$pseudotime.estimation$initiator.sample <- colnames(env$metacell.data)[1]
  }
  
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
  adj.matrix.kNNG[which(adj.matrix.kNNG<0)] <- 0
  
  if( any( is.infinite( shortest.paths(graph.adjacency(adj.matrix.kNNG,weighted=TRUE),colnames(env$metadata)[1]) ) ) )
  {
    util.warn("Disabling pseudotime estimation due to isolated subgraphs. Increase k.")
    env$preferences$pseudotime.estimation <- NULL
    return(env)
  }
  
  progressbar <- newProgressBar(min = 0, max = env$preferences$pseudotime.estimation$n.iterations); cat("\r")
  
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
    }
    
    T.list[[ng.i]] <- traj.it
    
    setTxtProgressBar( progressbar, progressbar$getVal()+1 )
  }
  
  progressbar$kill()
  
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
  env$pseudotime.trajectory <- (env$pseudotime.trajectory-min(env$pseudotime.trajectory))/(max(env$pseudotime.trajectory)-min(env$pseudotime.trajectory))
  return(env)
}


