pipeline.genesetStatisticIntegral <- function(env)
{
  spot.fisher.p <- function(spot)
  {
    spot$Fisher.p <- GeneSet.Fisher(unique(env$gene.info$ensembl.mapping$ensembl_gene_id[ which(env$gene.info$ensembl.mapping[,1]%in%spot$genes) ]),
                                    unique(env$gene.info$ensembl.mapping$ensembl_gene_id), env$gs.def.list, sort=TRUE)
    
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
