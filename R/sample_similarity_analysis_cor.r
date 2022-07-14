pipeline.sampleSimilarityAnalysisCor <- function(env)
{
  mixColors<-function(c1, c2,alpha=122 )
  {
    do.call( rgb,
             as.list( c( colMeans( sapply(c(2,4,6),function(i)
               c( strtoi(substr(c1,i,i+1),16), strtoi(substr(c2,i,i+1),16) )
               
             ) ), alpha ) /255 )
    )
  }
  
  pcm.module <- cor( env[[paste("spot.list.",env$preferences$standard.spot.modules,sep="")]]$spotdata )
  pcm.metadata <- cor( env$metadata )
  
  #### Correlation Spanning Tree ####
  
  filename <- file.path(paste(env$files.name, "- Results"), "Sample Similarity Analysis", "Correlation Spanning Tree.pdf")
  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 21/2.54, useDingbats=FALSE)
  
  for (i in 1:2 )
  {
    adj.matrix <- if( i == 1 ) pcm.module else pcm.metadata
    adj.matrix[ which(is.na(adj.matrix)) ] <- 0
    n <- if( i == 1 ) paste( "module data (", env$preferences$standard.spot.modules, ")") else "metagene data"
    
    g <- graph.adjacency(-adj.matrix, weighted=TRUE, mode="undirected")
    stg <- minimum.spanning.tree(g)
    E(stg)$weight <- 1
    layout <- layout_with_kk(stg)
    
    par(mar=c(1,1,1,1))
    plot(stg, layout=layout, vertex.size=5, vertex.label = rep("",ncol(adj.matrix)),
         vertex.color=env$group.colors, main=paste("Correlation Spanning Tree on",n))
      legend("bottomright", as.character(unique(env$group.labels)), cex=0.5, text.col=env$groupwise.group.colors, bg="white")
      box()
 
  }
  dev.off()

  
  #### Correlation Network ####
  
  filename <- file.path(paste(env$files.name, "- Results"), "Sample Similarity Analysis", "Correlation Network.pdf")
  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 21/2.54, useDingbats=FALSE)
  
  for (i in 1:2 )
  {
    adj.matrix <- if( i == 1 ) pcm.module else pcm.metadata
    adj.matrix[ which(is.na(adj.matrix)) ] <- 0
    n <- if( i == 1 ) paste( "module data (", env$preferences$standard.spot.modules, ")") else "metagene data"
    
    diag(adj.matrix) <- 0
    adj.matrix[which(adj.matrix < 0.5)] <- 0
    
    if (max(adj.matrix) > 0)
    {
      g <- graph.adjacency(adj.matrix, weighted=TRUE,  mode="undirected")
      E(g)$weight <- (2 + E(g)$weight)/2
      layout <- layout_with_kk( g )
      
      par(mar=c(1,1,1,1))      
      E(g)$color <- apply( get.edgelist( g ), 1, function(x) mixColors( env$group.colors[x[1]], env$group.colors[x[2]], alpha=40 ) )
      plot(g, layout=layout, vertex.size=ifelse(ncol(adj.matrix)<250, 5, 3),
           vertex.label = rep("",ncol(adj.matrix)),
           vertex.color=env$group.colors, main=paste("Correlation Network on",n))
      legend("bottomright", as.character(unique(env$group.labels)), cex=0.5, text.col=env$groupwise.group.colors, bg="white")
      box()
      
      if (ncol(adj.matrix) < 1000)
      {
        plot(g, layout=layout, vertex.size=ifelse(ncol(adj.matrix) < 250, 5, 3),
             vertex.label=colnames(adj.matrix),
             vertex.label.cex=if (ncol(adj.matrix)<100) 1.2 else 0.6,
             vertex.color=env$group.colors, main=paste("Correlation Network on",n))
        legend("bottomright", as.character(unique(env$group.labels)), cex=0.5, text.col=env$groupwise.group.colors, bg="white")
        box()
      }
    }
  }
  dev.off()
  
}
