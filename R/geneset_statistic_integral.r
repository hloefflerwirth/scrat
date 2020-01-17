pipeline.genesetStatisticIntegral <- function()
{
  spot.fisher.p <- function(spot)
  {
    spot$Fisher.p <- GeneSet.Fisher(unique(na.omit(gene.info$ids[spot$genes])),
                                    unique.protein.ids, gs.def.list, sort=TRUE)

    return(spot)
  }

  spot.list.overexpression$spots <<- lapply( spot.list.overexpression$spots, spot.fisher.p)
  spot.list.kmeans$spots <<- lapply( spot.list.kmeans$spots, spot.fisher.p)
  if (length(unique(group.labels)) > 1)
  {
    spot.list.group.overexpression$spots <<- lapply( spot.list.group.overexpression$spots, spot.fisher.p)
  }  
  spot.list.dmap$spots <<- lapply( spot.list.dmap$spots, spot.fisher.p)

}
