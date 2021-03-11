scrat.markerReport <- function(env,f)
{
  if(missing(env)||missing(f))
  {
    util.fatal("No environment or file supplied!")
    return(FALSE)
  }

  marker.set <- unlist(read.table(f,stringsAsFactors=FALSE))
  marker.set <- marker.set[which(marker.set!="")]
  n <- sub(".txt","",tail(strsplit(f,"/")[[1]],1))
  
  ylim <- quantile(env$seuratObject@assays$RNA@data,c(0.01,0.99))

  dir.create(paste(env$files.name, "- Results/Marker Set Reports"), showWarnings=FALSE)
  filename <- file.path(paste(env$files.name, "- Results"),"Marker Set Reports",paste(n,".pdf",sep=""))
  pdf(filename, 29.7/2.54, 21/2.54)
    
    #### Profile + Heatmap
    layout(matrix(c(1,2,3),ncol=1,byrow=T),heights=c(1.5,0.5,4))
    
    gs.indata <- env$seuratObject@assays$RNA@data[marker.set,]
    
    o.genes <- if(length(marker.set)>1) hclust(dist(gs.indata))$order else 1
    o.samples <- order(colMeans(gs.indata))
    
    
    offset <- min(colMeans(gs.indata))
    
    par(mar=c(0,10,5,16))
    barplot( colMeans(gs.indata)[o.samples]-offset, ylab="", xlab="", ylim=ylim-offset, xaxt="n", xaxs="i",
             col=env$group.colors[o.samples], xpd=TRUE, space=c(0,0), offset=0, axes=FALSE, border=if (ncol(env$seuratObject) < 100) "black" else NA )
    abline( h=0-offset, lty=2, col="gray" )
    axis( 2, at=seq(from=0,to=(max(colMeans(gs.indata))-offset),length.out=4), labels=round(seq(from=offset,to=max(colMeans(gs.indata)),length.out=4),2),las=2, cex.axis=1.4)
    mtext( n, side=3, cex=1.5, line=1 )
    mtext( "< e >", side=2, line=4.5, cex=1.25)
    
    
    par(new=T,mar=c(0,0,0,0))
    frame()
    legend(x=0.86,y=0.7,names(env$groupwise.group.colors),text.col=env$groupwise.group.colors)
    
    par(mar=c(0.5,10,0.5,16))
    image( matrix(1:ncol(env$seuratObject),ncol(env$seuratObject),1), col=env$group.colors[o.samples], axes=FALSE )
    box()
    
    par(mar=c(10,10,0,16))
    image(x=1:ncol(gs.indata), y=1:nrow(gs.indata),z=t(gs.indata[o.genes,o.samples,drop=FALSE]),col=env$color.palette.heatmaps(1000), axes=FALSE,zlim=max(max(gs.indata),-min(gs.indata))*c(-1,1),xlab="", ylab="")
    box()
    axis(2, 1:nrow(gs.indata), labels=rownames(gs.indata)[o.genes], las=2, line=-0.5, tick=0, cex.axis=0.9 )
    
    par(new=T,mar=c(40,74,0,1))
    image(matrix(c(1:1000), 1000, 1), axes=FALSE, col=env$color.palette.heatmaps(1000))
    box()
    axis( 1, at=c(0,1), labels=round(range(gs.indata),1), cex.axis=1.4 )
    mtext( bquote(Delta ~ "e"), side=1, cex=1.25, line=1 )  
    
    
    
    #### Profiles + Mapping
    layout(matrix(c(1,0,2,0),ncol=2),widths=c(1,0.5), heights=c(1,0.08))
    
    # barplot
    par(mar=c(14,4,4,1))
    barplot(colMeans(gs.indata), beside=TRUE,
            las=2, cex.names=1.2, cex.axis=1.4, col=env$group.colors, 
            ylim=ylim, border=if (ncol(env$seuratObject) < 100) "black" else NA,
            names.arg=rep("",ncol(env$seuratObject)))
    
    mtext( n, side=3, cex=1.5, line=1 )
    mtext( "< e >", side=2, line=2.5, cex=1.25)
    
    
    # population map
    par(mar=c(16, 0, 6, 3.5))
    
    n.map <- matrix(0,env$preferences$dim.1stLvlSom,env$preferences$dim.1stLvlSom)
    gs.nodes <- env$som.result$nodes[marker.set]
    n.map[as.numeric(names(table(gs.nodes)))] <- table(gs.nodes)
    n.map[which(n.map==0)] <- NA
    n.map <- matrix(n.map, env$preferences$dim.1stLvlSom)
    
    lim <- c(1,env$preferences$dim.1stLvlSom) + env$preferences$dim.1stLvlSom * 0.01 * c(-1, 1)
    colr <- env$color.palette.heatmaps(1000)[(na.omit(as.vector(n.map)) - min(n.map,na.rm=TRUE)) /
                                           max(1, (max(n.map,na.rm=TRUE) - min(n.map,na.rm=TRUE))) *
                                           999 + 1]
    
    plot(which(!is.na(n.map), arr.ind=TRUE), xlim=lim, ylim=lim, pch=16, axes=FALSE,
         xlab="",ylab="", xaxs="i", yaxs="i", col=colr,
         cex=0.5 + na.omit(as.vector(n.map)) / max(n.map,na.rm=TRUE) * 2.8)
    
    title(sub=paste("# features =", length(marker.set)),line=0)
    box()
    
    par(new=TRUE, mar=c(32,21.5,6,0.5))
    image(matrix(1:100, 1, 100), col = env$color.palette.heatmaps(1000), axes=FALSE)
    axis(2, at=c(0,1), c(min(n.map,na.rm=TRUE), max(n.map,na.rm=TRUE)), las=2, tick=FALSE, pos=-0.36)
    box()
    
    
    #### Profiles - Group boxplots
    layout(1)
    par(mar=c(12,6,4,5))
    
    mean.boxes <- by(colMeans(gs.indata), env$group.labels, c)[unique(env$group.labels)]
    boxplot(mean.boxes, col=env$groupwise.group.colors, ylim=ylim+sum(abs(ylim))*0.1*c(-1,1), axes=FALSE, yaxs="i")
    box()
    axis(1, seq_along(unique(env$group.labels)), unique(env$group.labels), las=2, tick=FALSE)
    axis(2, las=2)
    
    abline(h=0, lty=2)
    
    mtext( n, side=3, cex=1.5, line=1 )
    mtext("< e >", side=2, line=2.5, cex=1.25)
    
  dev.off()
  
  return(invisible(TRUE))
}
