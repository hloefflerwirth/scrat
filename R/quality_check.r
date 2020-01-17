pipeline.qualityCheck <- function()
{
  dir.create(paste(files.name, "- Results/Data Overview"), showWarnings=FALSE)

  mi = min(indata); ma = max(indata)
  densities = apply(indata,2,function(x) density(x, from=mi, to=ma))
  densities.y = t(sapply(densities, function(x) x$y))
  densities.x = densities[[1]]$x

  filename <- file.path(paste(files.name, "- Results"), "Data Overview", "Data Distribution.pdf")
  util.info("Writing:", filename)
  pdf(filename, 29.7/2.54, 21/2.54, useDingbats=FALSE)

  plot(densities.x, densities.y[1,], main="Input data distribution", xlim=range(indata), ylim=range(densities.y), type="l", xlab="log Expression", ylab="Density")
    dummy=sapply(densities, lines)

	
	indata.sample.sd <- apply(indata, 2, sd)
	indata.sample.mean <- colMeans(indata)

	Q13 <- quantile( indata.sample.mean, c(0.25,0.75) )
	IQR1.mean <- c( Q13[1] - 1*diff( Q13 ), Q13[2] + 1*diff( Q13 ) )
	IQR3.mean <- c( Q13[1] - 3*diff( Q13 ), Q13[2] + 3*diff( Q13 ) )
	Q13 <- quantile( indata.sample.sd, c(0.25,0.75) )
	IQR1.sd <- c( Q13[1] - 1*diff( Q13 ), Q13[2] + 1*diff( Q13 ) )
	IQR3.sd <- c( Q13[1] - 3*diff( Q13 ), Q13[2] + 3*diff( Q13 ) )

	outlier <- names( which( indata.sample.mean < IQR1.mean[1] |
														 indata.sample.mean > IQR1.mean[2] |
														 indata.sample.sd < IQR1.sd[1] |
														 indata.sample.sd > IQR1.sd[2] ) )

	par(mfrow=c(1,1), mar=c(5,4,3,2))
	plot( indata.sample.mean, indata.sample.sd, pch=16, col=group.colors, xlab="mean expression level", ylab="standard deviation of expression" )
		abline( v=IQR1.mean, col="gray20", lty=2 )
		abline( v=IQR3.mean, col="gray20", lty=3 )
		abline( h=IQR1.sd, col="gray20", lty=2 )
		abline( h=IQR3.sd, col="gray20", lty=3 )
		legend("topright",c("1x IQR","3x IQR"),lty=c(2,3),col="gray20")
	if(length(outlier)>0) text(indata.sample.mean[outlier],indata.sample.sd[outlier]+diff(range(indata.sample.sd))*0.017,outlier)
		
		
  par(mfrow=c(2,1),mar=c(5,3,3,2))

  barplot(indata.sample.mean, col=group.colors, main="Sample mean expression", names.arg=if (ncol(indata)<80) colnames(indata) else rep("",ncol(indata)), las=2, cex.main=1, cex.lab=1, cex.axis=0.8, cex.names=0.6, border = ifelse(ncol(indata)<80,"black",NA), ylim=range(indata.sample.mean)*c(0.99,1.01), xpd=FALSE)
    box()
	

  if (length(unique(group.labels)) > 1)
  {
    mean.boxes = by(indata.sample.mean, group.labels, c)[unique(group.labels)]
    par(mar=c(5,3,0,2))
    boxplot(mean.boxes, col=groupwise.group.colors, las=2, main="", cex.main=1, cex.axis=0.8, xaxt="n")
    axis(1, seq_along(groupwise.group.colors), unique(group.labels), las=2, cex.axis=0.8)
  }


  dev.off()
}
