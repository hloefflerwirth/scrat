pipeline.sampleExpressionPortraits <- function(env)
{
  ## Expression portraits
  filename <- file.path(paste(env$files.name, "- Results"), "Expression Portraits.pdf")
  util.info("Writing:", filename)

  pdf(filename, 29.7/2.54, 21/2.54, useDingbats=FALSE)
  par(mfrow=c(7, 12))
  par(mar=c(0.3, 0.9, 4.5, 0.9))
  count.col <- 0

  for (gl in seq_along(unique(env$group.labels)))
  {
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1))

    if (length(unique(env$group.labels)) > 1)
    {
      mtext(unique(env$group.labels)[gl], side=3, line=2, cex=1.5, at=0, font=3,
            adj=0, col=env$groupwise.group.colors[gl])
    }

    par(new=TRUE)

    for (j in which(env$group.labels == unique(env$group.labels)[gl]))
    {
      image(matrix(env$metadata[,j], env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom),
            axes=FALSE, col=env$color.palette.portraits(1000))

      title(paste(j,":",colnames(env$indata)[j]), line=1, cex.main=0.8)

      if (length(unique(env$group.labels)) > 1)
      {
        title(paste("S =",round(env$group.silhouette.coef[j],2)), line=0.2,
              cex.main=0.6, col.main=env$groupwise.group.colors[gl])
      }

      box()
      count.col <- count.col + 1
    }

    if (count.col %% 12 != 0)
    {
      for (j in 1:(12 - count.col %% 12))
      {
        plot(0, type="n", axes=FALSE, xlab="", ylab="")
        count.col <- count.col + 1
      }
    }
  }

  dev.off()

}
