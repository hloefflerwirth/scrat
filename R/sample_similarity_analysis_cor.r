pipeline.sampleSimilarityAnalysisCor <- function(env)
{
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

}
