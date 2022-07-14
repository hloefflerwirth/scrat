modules.CSV.sheets <- function(env, spot.list, main, path)
{
  
  for (m in seq_along(spot.list$spots))
  {
    basename <- paste(main, " ", names(spot.list$spots)[m], ".csv", sep="")

    if (length(spot.list$spots[[m]]$genes) <= 0 || length(spot.list$spots[[m]]$genes) > 2000)
    {
      next
    }

    ## CSV Table

    o <- spot.list$spots[[m]]$genes
    
    out <- data.frame(Rank=c(seq_along(spot.list$spots[[m]]$genes)),
                      ID=o,
                      Symbol=env$gene.info$names[o])


    out <- cbind(out,
                    "Metagene"=env$gene.info$coordinates[o],
                    "Chromosome"=paste( env$gene.info$chr.name[o], env$gene.info$chr.band[o]),
                    "Description"=env$gene.info$descriptions[o])

    env$csv.function(out, file.path(path, basename))
  }
}
