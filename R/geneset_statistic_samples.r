pipeline.genesetStatisticSamples <- function(env)
{

  ### perform GS analysis ###
  env$t.ensID.m <- env$t.g.m[which(rownames(env$indata) %in% names(env$gene.info$ids)),]
  env$t.ensID.m <- do.call(rbind, by(env$t.ensID.m, env$gene.info$ids, colMeans))

  if (env$preferences$activated.modules$geneset.analysis.exact)
  {
    gs.null.list <- list()

    for (i in seq_along(env$gs.def.list))
    {
      gs.null.list[[i]] <-
        list(Genes=sample(env$unique.protein.ids, length(env$gs.def.list[[i]]$Genes)))
    }

    null.scores <- sapply( 1:ncol(env$indata), function(m)
    {
      all.gene.statistic <- env$t.ensID.m[,m]
      spot.gene.ids <- env$unique.protein.ids

      scores <- GeneSet.GSZ(spot.gene.ids, all.gene.statistic, gs.null.list)

      return(scores)
    })

    null.culdensity <- ecdf(abs(unlist(null.scores)))
  }

  env$spot.list.samples <- lapply(seq_along(env$spot.list.samples) , function(m)
  {
    x <- env$spot.list.samples[[m]]
    all.gene.statistic <- env$t.ensID.m[,m]
    spot.gene.ids <- env$unique.protein.ids

    x$GSZ.score <-
      GeneSet.GSZ(spot.gene.ids, all.gene.statistic, env$gs.def.list, sort=FALSE)
#      GeneSet.maxmean(all.gene.statistic, env$gs.def.list)

    if (env$preferences$activated.modules$geneset.analysis.exact)
    {
       x$GSZ.p.value <- 1 - null.culdensity(abs(x$GSZ.score))
       names(x$GSZ.p.value) <- names(x$GSZ.score)
    }

    return(x)
  })
  names(env$spot.list.samples) <- colnames(env$indata)
  
  ### GSZ table output ###
  env$samples.GSZ.scores <- do.call(cbind, lapply(env$spot.list.samples, function(x)
  {
    return(x$GSZ.score[names(env$gs.def.list)])
  }))
  return(env)
}
