pipeline.summarySheetsModules <- function()
{
  environment(modules.CSV.sheets) <- environment()
  environment(modules.report.sheets) <- environment()
  
  dirname <- paste(files.name, "- Results/Summary Sheets - Modules")
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  #### Overexpression Spots ####
  
  dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","Overexpression Spots" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)

  modules.report.sheets(spot.list=spot.list.overexpression, main="Overexpression Spots", path=file.path(dirname,"Report.pdf") )
  
  #### K-Means Cluster ####
  
  dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","K-Means Cluster" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  
  modules.report.sheets(spot.list=spot.list.kmeans, main="K-Means Cluster", path=file.path(dirname,"Report.pdf") )
  
  #### D-Clusters ####
  
  dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","D-Cluster" )
  util.info("Writing:", file.path(dirname, "*.pdf"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  
  modules.report.sheets(spot.list=spot.list.dmap, main="D-Cluster", path=file.path(dirname,"Report.pdf") )
  
  #### Group Overexpression Spots ####
  
  if (length(unique(group.labels)) > 1)
  {
    dirname <- file.path( paste(files.name, "- Results"),"Summary Sheets - Modules","Group Overexpression Spots" )
    util.info("Writing:", file.path(dirname, "*.pdf"))
    dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
    
    modules.report.sheets(spot.list=spot.list.group.overexpression, main="Group Overexpression Spots", path=file.path(dirname,"Report.pdf") )
  }

  #### module gene lists CSV sheets ####
   
  dirname <- file.path(output.paths["CSV"], "Spot Lists")
  util.info("Writing:", file.path(dirname, "*.csv"))
  dir.create(dirname, showWarnings=FALSE, recursive=TRUE)
  
  modules.CSV.sheets(spot.list=spot.list.overexpression, main="Overexpression Spots", path=dirname )
  modules.CSV.sheets(spot.list=spot.list.kmeans, main="K-Means Cluster", path=dirname )
  modules.CSV.sheets(spot.list=spot.list.dmap, main="D-Cluster", path=dirname )
  
  if (length(unique(group.labels)) > 1)
  {
    modules.CSV.sheets(spot.list=spot.list.group.overexpression, main="Group Overexpression Spots", path=dirname )
  }
  
}