pipeline.calcStatistics <- function(env)
{
  util.info("Calculating Single Gene Statistic")
  progressbar <- newProgressBar(min = 0, max = ncol(env$indata)); cat("\r")
 
  # Calculate T-score and significance
  
  sd.g.m <- matrix(NA, nrow(env$indata), ncol(env$indata), dimnames=list(rownames(env$indata), colnames(env$indata)))
  
  t.g.m <- matrix(NA, nrow(env$indata), ncol(env$indata), dimnames=list(rownames(env$indata), colnames(env$indata)))
  env$p.g.m <- matrix(NA, nrow(env$indata), ncol(env$indata), dimnames=list(rownames(env$indata), colnames(env$indata)))
  
  n.0.m <- rep(NA, ncol(env$indata))
  names(n.0.m) <- colnames(env$indata)
  
  env$perc.DE.m <- rep(NA, ncol(env$indata))
  names(env$perc.DE.m) <- colnames(env$indata)
  
  env$fdr.g.m <- matrix(NA, nrow(env$indata), ncol(env$indata), dimnames=list(rownames(env$indata), colnames(env$indata)))
  
  
  for (m in 1:ncol(env$indata))
  {
    env$p.g.m[,m] <- apply( env$indata, 1, function(x)
    {
      if( all(x[-m] == x[-m][1]) ) return(1) 
      
      return( t.test( x[m], x[-m], var.equal=TRUE )$p.value )
    } )
    
    suppressWarnings({ try.res <- try({
      fdrtool.result <- fdrtool(env$p.g.m[,m], statistic="pvalue", verbose=FALSE, plot=FALSE)
    }, silent=TRUE) })
    
    if (!is(try.res,"try-error"))
    {
      env$fdr.g.m[,m] <- fdrtool.result$lfdr
      
      n.0.m[m] <- fdrtool.result$param[1,"eta0"]
      env$perc.DE.m[m] <- 1 - n.0.m[m]
    } else # happens for eg phenotype data
    {
      env$p.g.m[,m] <- order(env$indata[,m]) / nrow(env$indata)
      env$fdr.g.m[,m] <- env$p.g.m[,m]
      
      env$perc.DE.m[m] <- 0
    }
    
    setTxtProgressBar( progressbar, progressbar$getVal()+1 )
  }
  
  progressbar$kill()

  return(env)
}
