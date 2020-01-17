pipeline.sampleSimilarityAnalysisCor <- function()
{
  pcm.module <- cor( get(paste("spot.list.",preferences$standard.spot.modules,sep=""))$spotdata )
  pcm.metadata <- cor( metadata )
  
  #### Correlation Spanning Tree ####
  
  filename <- file.path(paste(files.name, "- Results"), "Sample Similarity Analysis", "Correlation Spanning Tree.pdf")
  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 21/2.54, useDingbats=FALSE)
  
  for (i in 1:2 )
  {
    adj.matrix <- if( i == 1 ) pcm.module else pcm.metadata
    n <- if( i == 1 ) paste( "module data (", preferences$standard.spot.modules, ")") else "metagene data"
    
    g <- graph.adjacency(-adj.matrix, weighted=TRUE, mode="undirected")
    stg <- minimum.spanning.tree(g)
    E(stg)$weight <- 1
    layout <- layout_with_kk(stg)
    
    par(mar=c(1,1,1,1))
    plot(stg, layout=layout, vertex.size=5, vertex.label = rep("",ncol(adj.matrix)),
         vertex.color=group.colors, main=paste("Correlation Spanning Tree on",n))
      legend("bottomright", as.character(unique(group.labels)), cex=0.5, text.col=groupwise.group.colors, bg="white")
      box()
 
  }
  dev.off()

}
