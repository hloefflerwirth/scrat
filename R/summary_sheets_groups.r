pipeline.summarySheetsGroups <- function(env)
{
  group.metadata <- do.call(cbind, by(t(env$metadata), env$group.labels, colMeans))[,unique(env$group.labels)]
  group.metadata.sd <- do.call(cbind, by(t(env$metadata), env$group.labels, function(x) { apply(x,2,sd) }))[,unique(env$group.labels)]

  loglog.group.metadata <- apply(group.metadata, 2, function(x)
  {
    meta.sign <- sign(x)
    meta <- log10(abs(x))
    meta <- meta - min(meta, na.rm=TRUE)
    return(meta * meta.sign)
  })
  WAD.group.metadata <- apply(group.metadata ,2, function(x) { x * ((x - min(x)) / (max(x) - min(x))) })

  bleached.group.metadata <- group.metadata
  bleached.WAD.group.metadata <- WAD.group.metadata
  bleached.loglog.group.metadata <- loglog.group.metadata

  for (i in seq_along(unique(env$group.labels)))
  {
    pos.metagenes <- which(group.metadata[,i] >= 0)
    neg.metagenes <- which(group.metadata[,i] < 0)

    bleached.group.metadata[pos.metagenes,i] <-
      bleached.group.metadata[pos.metagenes,i] -
      pmin(bleached.group.metadata[pos.metagenes,i],
           apply(group.metadata[pos.metagenes,-i,drop=FALSE], 1, max))

    bleached.group.metadata[neg.metagenes,i] <-
      bleached.group.metadata[neg.metagenes,i] -
      pmax(bleached.group.metadata[neg.metagenes,i],
           apply(group.metadata[neg.metagenes,-i,drop=FALSE], 1, min))

    pos.metagenes <- which(WAD.group.metadata[,i] >= 0)
    neg.metagenes <- which(WAD.group.metadata[,i] < 0)

    bleached.WAD.group.metadata[pos.metagenes,i] <-
      bleached.WAD.group.metadata[pos.metagenes,i] -
      pmin(bleached.WAD.group.metadata[pos.metagenes,i],
           apply(WAD.group.metadata[pos.metagenes,-i,drop=FALSE], 1, max))

    bleached.WAD.group.metadata[neg.metagenes,i] <-
      bleached.WAD.group.metadata[neg.metagenes,i] -
      pmax(bleached.WAD.group.metadata[neg.metagenes,i],
           apply(WAD.group.metadata[neg.metagenes,-i,drop=FALSE], 1, min))

    pos.metagenes <- which(loglog.group.metadata[,i] >= 0)
    neg.metagenes <- which(loglog.group.metadata[,i] < 0)

    bleached.loglog.group.metadata[pos.metagenes,i] <-
      bleached.loglog.group.metadata[pos.metagenes,i] -
      pmin(bleached.loglog.group.metadata[pos.metagenes,i],
           apply(loglog.group.metadata[pos.metagenes,-i,drop=FALSE], 1, max))

    bleached.loglog.group.metadata[neg.metagenes,i] <-
      bleached.loglog.group.metadata[neg.metagenes,i] -
      pmax(bleached.loglog.group.metadata[neg.metagenes,i],
           apply(loglog.group.metadata[neg.metagenes,-i,drop=FALSE], 1, min))

  }


  ###### Group stability scores
  filename <- file.path(paste(env$files.name, "- Results"),
                        "Summary Sheets - Groups",
                        "Group Assignment.pdf")

  util.info("Writing:", filename)
  pdf(filename, 21/2.54, 29.7/2.54, useDingbats=FALSE)


  S <- tapply( env$group.silhouette.coef, env$group.labels, sort, decreasing=TRUE, simplify=FALSE )[unique(env$group.labels)]
  names(S) <- NA
  S <- unlist(S)
  names(S) <- sub( paste("^NA.",sep=""), "", names(S) )
  
  
  PCM <- cor( env$metadata )
  diag(PCM) <- NA
   
  group.correlations <- sapply( seq(ncol(env$metadata)), function(i)
  {
    mean.group.correlations <- tapply( PCM[,i], env$group.labels, mean, na.rm=TRUE )[unique(env$group.labels)]
    
    return(  mean.group.correlations )
  } )
  colnames(group.correlations) <- colnames(env$indata)
  group.correlations[which(is.nan(group.correlations))] <- 0
    
  
  #layout(matrix(c(0,1,2,0),ncol=1))
  par(mfrow=c(2,1))
  par(mar=c(5,3,3,2))
  
  b<-barplot( S, col=env$group.colors[names(S)], main="Correlation Silhouette", names.arg=if(ncol(env$indata)<80) names(S) else rep("",length(S)), las=2, cex.main=1, cex.lab=1, cex.axis=0.8, cex.names=0.6, border = ifelse(ncol(env$indata)<80,"black",NA), xpd=FALSE, ylim=c(-.25,1) )  
  mtext("S",2,line=1.9,cex=0.8)
  abline( h=c(0,0.25,0.5,0.75), lty=2, col="gray80" )
  title( main= bquote("<" ~ s ~ "> = " ~ .(round(mean(S),2))), line=0.5, cex.main=1 )
  box()
  points( b, rep(-0.2,ncol(env$indata)), pch=15, cex=1, col=env$groupwise.group.colors[apply( group.correlations[,names(S)], 2, which.max )] )
  
  mean.boxes <- by( S, env$group.labels, c )[ unique( env$group.labels ) ]
  mean.mean.S <- sapply( mean.boxes, mean )
  
  par(mar=c(5,3,0,2))
  boxplot( mean.boxes, col=env$groupwise.group.colors, las=2, main="", cex.main=1, cex.axis=0.8, xaxt="n", ylim=c(-.25,1) )
  mtext("S",2,line=1.9,cex=0.8)
  abline( h=c(0,0.25,0.5,0.75), lty=2, col="gray80" )
  axis( 1, 1:length(env$groupwise.group.colors), paste( unique(env$group.labels), "\n<s> =", round(mean.mean.S,2) ), las=2, cex.axis=0.8 )
  
  dev.off()

}
