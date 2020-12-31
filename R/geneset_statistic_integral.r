pipeline.genesetStatisticIntegral <- function(env)
{
  spot.fisher.p <- function(spot)
  {
    spot$Fisher.p <- GeneSet.Fisher(unique(na.omit(env$gene.info$ids[env$spot$genes])),
                                    env$unique.protein.ids, env$gs.def.list, sort=TRUE)

    return(spot)
  }

  env$spot.list.overexpression$spots <- lapply( env$spot.list.overexpression$spots, spot.fisher.p)
  env$spot.list.kmeans$spots <- lapply( env$spot.list.kmeans$spots, spot.fisher.p)
  if (length(unique(env$group.labels)) > 1)
  {
    env$spot.list.group.overexpression$spots <- lapply( env$spot.list.group.overexpression$spots, spot.fisher.p)
  }  
  env$spot.list.dmap$spots <- lapply( env$spot.list.dmap$spots, spot.fisher.p)

  return(env)
}
