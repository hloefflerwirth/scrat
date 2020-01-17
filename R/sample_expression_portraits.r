pipeline.sampleExpressionPortraits <- function()
{
  ## Expression portraits
  filename <- file.path(paste(files.name, "- Results"), "Expression Portraits.pdf")
  util.info("Writing:", filename)

  pdf(filename, 29.7/2.54, 21/2.54, useDingbats=FALSE)
  par(mfrow=c(7, 12))
  par(mar=c(0.3, 0.9, 4.5, 0.9))
  count.col <- 0

  for (gl in seq_along(unique(group.labels)))
  {
    plot(0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(0,1))

    if (length(unique(group.labels)) > 1)
    {
      mtext(unique(group.labels)[gl], side=3, line=2, cex=1.5, at=0, font=3,
            adj=0, col=groupwise.group.colors[gl])
    }

    par(new=TRUE)

    for (j in which(group.labels == unique(group.labels)[gl]))
    {
      image(matrix(metadata[,j], preferences$dim.1stLvlSom, preferences$dim.1stLvlSom),
            axes=FALSE, col=color.palette.portraits(1000))

      title(paste(j,":",colnames(indata)[j]), line=1, cex.main=0.8)

      if (length(unique(group.labels)) > 1)
      {
        title(paste("S =",round(group.silhouette.coef[j],2)), line=0.2,
              cex.main=0.6, col.main=groupwise.group.colors[gl])
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
