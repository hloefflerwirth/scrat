pipeline.checkInputParameters <- function(env)
{
  #### check preferences ####
  if (!is.character(env$preferences$dataset.name))
  {
    util.warn("Invalid value of \"dataset.name\". Using \"Unnamed\"")
    env$preferences$dataset.name <- "Unnamed"
  }

  if ( env$preferences$dim.1stLvlSom!="auto" && !is.numeric(env$preferences$dim.1stLvlSom) || env$preferences$dim.1stLvlSom < 1)
  {
    util.warn("Invalid value of \"dim.1stLvlSom\". Using size recommendation")
    env$preferences$dim.1stLvlSom <- "auto"
  }

  if (!is.numeric(env$preferences$training.extension) ||
      env$preferences$training.extension < 1 ||
      env$preferences$training.extension > 10)
  {
    util.warn("Invalid value of \"training.extension\". Using 1")
    env$preferences$training.extension <- 1
  }

  if (!is.numeric(env$preferences$rotate.SOM.portraits) ||
      env$preferences$rotate.SOM.portraits < 0 ||
      env$preferences$rotate.SOM.portraits > 4)
  {
    util.warn("Invalid value of \"rotate.SOM.portraits\". Using 0")
    env$preferences$rotate.SOM.portraits <- 0
  }

  if (!is.logical(env$preferences$flip.SOM.portraits))
  {
    util.warn("Invalid value of \"flip.SOM.portraits\". Using FALSE")
    env$preferences$flip.SOM.portraits <- FALSE
  }

  if (!is.character(env$preferences$database.biomart))
  {
    util.warn("Invalid value of \"database.biomart\". Using \"\"")
    env$preferences$database.biomart <- ""
  }

  if (!is.character(env$preferences$database.host))
  {
    util.warn("Invalid value of \"database.host\". Using \"\"")
    env$preferences$database.host <- ""
  }

  if (!is.character(env$preferences$database.dataset))
  {
    util.warn("Invalid value of \"database.dataset\". Using \"\"")
    env$preferences$database.dataset <- ""
  }

  if (!is.character(env$preferences$database.id.type))
  {
    util.warn("Invalid value of \"database.id.type\". Using \"\"")
    env$preferences$database.id.type <- ""
  }

  if (!is.list(env$preferences$activated.modules))
  {
    util.warn("Invalid value of \"activated.modules\". Using all analysis modules")
    env$preferences$activated.modules = list( "reporting" = TRUE,
                              "primary.analysis" = TRUE,
                              "sample.similarity.analysis" = TRUE,
                              "geneset.analysis" = TRUE,
                              "group.analysis" = TRUE,
                              "difference.analysis" = TRUE )
  } else
  {
    if (!is.logical(env$preferences$activated.modules$reporting))
    {
      util.warn("Invalid value of \"activated.modules$reporting\". Using TRUE")
      env$preferences$activated.modules$reporting <- TRUE
    }
    if (!is.logical(env$preferences$activated.modules$primary.analysis))
    {
      util.warn("Invalid value of \"activated.modules$primary.analysis\". Using TRUE")
      env$preferences$activated.modules$primary.analysis <- TRUE
    } else
    if (!env$preferences$activated.modules$primary.analysis && is.null(env$som.result))
    {
      util.warn("No primary analysis perfomed yet. Setting \"activated.modules$primary.analysis\" to TRUE")
      env$preferences$activated.modules$primary.analysis <- TRUE
    }
    if (!is.logical(env$preferences$activated.modules$sample.similarity.analysis))
    {
      util.warn("Invalid value of \"activated.modules$sample.similarity.analysis\". Using TRUE")
      env$preferences$activated.modules$sample.similarity.analysis <- TRUE
    }
    if (!is.logical(env$preferences$activated.modules$geneset.analysis))
    {
      util.warn("Invalid value of \"activated.modules$geneset.analysis\". Using TRUE")
      env$preferences$activated.modules$geneset.analysis <- TRUE
    }
    if (!is.logical(env$preferences$activated.modules$group.analysis))
    {
      util.warn("Invalid value of \"activated.modules$group.analysis\". Using TRUE")
      env$preferences$activated.modules$group.analysis <- TRUE
    }
    if (!is.logical(env$preferences$activated.modules$difference.analysis))
    {
      util.warn("Invalid value of \"activated.modules$difference.analysis\". Using TRUE")
      env$preferences$activated.modules$difference.analysis <- TRUE
    }
  }

  if (!is.character(env$preferences$standard.spot.modules) || length(env$preferences$standard.spot.modules)!=1 ||
      !env$preferences$standard.spot.modules %in% c("overexpression","underexpression","kmeans","correlation","group.overexpression","dmap") )
  {
    util.warn("Invalid value of \"standard.spot.modules\". Using \"dmap\"")
    env$preferences$standard.spot.modules <- "dmap"
  }

  if (!is.numeric(env$preferences$spot.coresize.modules) ||
      env$preferences$spot.coresize.modules < 1 ||
      env$preferences$spot.coresize.modules > 20)
  {
    util.warn("Invalid value of \"spot.coresize.modules\". Using 3")
    env$preferences$spot.coresize.modules <- 3
  }

  if (!is.numeric(env$preferences$spot.threshold.modules) ||
      env$preferences$spot.threshold.modules <= 0 ||
      env$preferences$spot.threshold.modules >= 1)
  {
    util.warn("Invalid value of \"spot.threshold.modules\". Using 0.95")
    env$preferences$spot.threshold.modules <- 0.95
  }

  if (!is.numeric(env$preferences$spot.coresize.groupmap) ||
      env$preferences$spot.coresize.groupmap < 1 ||
      env$preferences$spot.coresize.groupmap > 20)
  {
    util.warn("Invalid value of \"spot.coresize.groupmap\". Using 5")
    env$preferences$spot.coresize.groupmap <- 5
  }

  if (!is.numeric(env$preferences$spot.threshold.groupmap) ||
      env$preferences$spot.threshold.groupmap <= 0 ||
      env$preferences$spot.threshold.groupmap >= 1)
  {
    util.warn("Invalid value of \"spot.threshold.groupmap\". Using 0.75")
    env$preferences$spot.threshold.groupmap <- 0.75
  }

  if (!is.character(env$preferences$dim.reduction) || length(env$preferences$dim.reduction)!=1 ||
      !env$preferences$dim.reduction %in% c("tsne", "pca", "umap") )
  {
    util.warn("Invalid value of \"dim.reduction\". Using \"tsne\"")
    env$preferences$dim.reduction <- "tsne"
  }

  if (!is.list(env$preferences$preprocessing))
  {
    util.warn("Invalid value of \"preprocessing\". Using default setting.")
    env$preferences$preprocessing = list(
      # count.processing = FALSE,
      cellcycle.correction = FALSE,
      create.meta.cell = FALSE,
      feature.centralization = TRUE,
      sample.quantile.normalization = TRUE )
      # seurat.normalize = TRUE,

  } else
  {
    if (!is.logical(env$preferences$preprocessing$cellcycle.correction))
    {
      util.warn("Invalid value of \"preprocessing$cellcycle.correction\". Using FALSE")
      env$preferences$preprocessing$cellcycle.correction <- FALSE
    }
    if (!is.logical(env$preferences$preprocessing$create.meta.cell))
    {
      util.warn("Invalid value of \"preprocessing$create.meta.cell\". Using FALSE")
      env$preferences$preprocessing$create.meta.cell <- FALSE
    }
    if (!is.logical(env$preferences$preprocessing$feature.centralization))
    {
      util.warn("Invalid value of \"preprocessing$feature.centralization\". Using TRUE")
      env$preferences$preprocessing$feature.centralization <- TRUE
    }
    if (!is.logical(env$preferences$preprocessing$sample.quantile.normalization))
    {
      util.warn("Invalid value of \"preprocessing$sample.quantile.normalization\". Using TRUE")
      env$preferences$preprocessing$sample.quantile.normalization <- TRUE
    }
  }

  #### check input data ####

  if(is.null(env$indata)){
    util.fatal("No data supplied!")
    env$passedInputChecking <- FALSE
    return(env)
  }

  if (is(env$indata, "character")) {
    util.info("Loading data from 10X from directory. This may take several time until next notification.")

    suppressWarnings({
      try.res <- try({
        env$indata <- Read10X(data.dir = env$indata)
        env$seuratObject <- CreateSeuratObject(counts = env$indata, project = env$preferences$dataset.name)
        env$indata <- NULL
      }, silent=TRUE)
    })

    if (is(try.res,"try-error")){
      util.fatal("Given diretory is invalid, cannot read input data.")
      env$passedInputChecking <- FALSE
      return(env)
    }

    # env$preferences$indata.counts <- TRUE
  } else
  if( is(env$indata,"Seurat") )
  {
    env$seuratObject <- env$indata
    env$indata <- NULL

    if( !"assays" %in% slotNames(env$seuratObject) ||
        !"RNA" %in% names(env$seuratObject@assays) ||
        !"counts" %in% slotNames(env$seuratObject@assays$RNA) )
    {
      env$seuratObject <- NULL
    }

  } else
  if( is(env$indata,"dgCMatrix") )
  {
    env$seuratObject <- suppressWarnings({ CreateSeuratObject(counts = env$indata, project = env$preferences$dataset.name) })
    env$indata <- NULL

  } else
  if( is(env$indata,"matrix") )
  {
    if( is.null(colnames(env$indata)) )
    {
      colnames(env$indata) <- paste( "cell", as.character(1:ncol(env$indata)) )
    }
    if( is.null(rownames(env$indata)) )
    {
      rownames(env$indata) <- paste( "gene", as.character(1:nrow(env$indata)) )
    }

    env$seuratObject <- suppressWarnings({ CreateSeuratObject(counts = env$indata, project = env$preferences$dataset.name) })
    env$indata <- NULL

  }

  if( !is(env$seuratObject,"Seurat") )
  {
    util.fatal("Data supplied in unknown format! Provide a Seurat v4 object, dgCMatrix, or simple matrix.")
    env$passedInputChecking <- FALSE
    return(env)
  }


  if ( "group.labels" %in% names( env$seuratObject ) )
  {
    env$group.labels <- env$seuratObject$group.labels

    if( "group.colors" %in% names( env$seuratObject ) )
      env$group.colors <- env$seuratObject$group.colors
  }

  if( length(env$group.labels)==1 && env$group.labels=="auto" )
  {
    env$group.labels <- rep("auto",ncol(env$seuratObject))
    names(env$group.labels) <- colnames(env$seuratObject)
  }

  # check group.labels and group.colors
  if ((!is.null(env$group.labels) && length(env$group.labels) != ncol(env$seuratObject)) ||
      (!is.null(env$group.colors) && length(env$group.colors) != ncol(env$seuratObject)))
  {
    env$group.labels <- NULL
    env$group.colors <- NULL
    util.warn("Group assignment doesnt fit number of samples")
  }
  if (!is.null(env$group.labels) && max(table(env$group.labels)) == 1)
  {
    env$group.labels <- NULL
    env$group.colors <- NULL
    util.warn("Each sample has an own group")
  }

  if (!is.null(env$group.labels))
  {
    for (sample in unique(colnames(env$seuratObject)))
    {
      if (length(unique(env$group.labels[which(colnames(env$seuratObject) == sample)])) > 1)
      {
        util.warn("Sample is in multiple groups:", sample)
        env$group.labels <- NULL
        env$group.colors <- NULL
        break
      }
    }
  }

  if (!is.null(env$group.labels))
  {
    env$group.labels <- as.character(env$group.labels)
    names(env$group.labels) <- colnames(env$seuratObject)

    if (is.null(env$group.colors))
    {
      env$group.colors <- rep("", ncol(env$seuratObject))

      for (i in seq_along(unique(env$group.labels)))
      {
        env$group.colors[which(env$group.labels == unique(env$group.labels)[i])] <-
          colorRampPalette(c("blue3", "blue", "green3", "gold", "red", "red3"))(length(unique(env$group.labels)))[i]
      }
    }

    # catch userdefined group.colors --> convert to #hex
    if (length(unique(substr(env$group.colors, 1, 1)) > 1) || unique(substr(env$group.colors, 1, 1))[1] != "#")
    {
      env$group.colors <- apply(col2rgb(env$group.colors), 2, function(x) { rgb(x[1]/255, x[2]/255, x[3]/255) })
    }
    names(env$group.colors) <- colnames(env$seuratObject)
  } else
  {
    # env$group.labels <- rep("auto",ncol(env$seuratObject))
    env$group.labels <- env$seuratObject$orig.ident
    names(env$group.labels) <- colnames(env$seuratObject)

    env$group.colors <- color.palette.discrete(length(unique(env$group.labels)))[match(env$group.labels, unique(env$group.labels))]
    names(env$group.colors) <- colnames(env$seuratObject)
  }

  env$groupwise.group.colors <- env$group.colors[match(unique(env$group.labels), env$group.labels)]
  names(env$groupwise.group.colors) <- unique(env$group.labels)

  env$seuratObject$group.labels <- env$group.labels
  env$seuratObject$group.colors <- env$group.colors



  if (env$preferences$dim.1stLvlSom == "auto")
  {
    n.sample.interval <- cut( ncol(env$seuratObject), breaks=c(0,100,500,1000,5000,Inf), labels=c(1:5) )
    n.feature.interval <- cut( nrow(env$seuratObject), breaks=c(0,1000,10000,Inf), labels=c(1:3) )
    recommendation <- matrix(c(seq(20,40,5),seq(30,50,5),seq(40,60,5)),nrow=3,byrow=TRUE)

    env$preferences$dim.1stLvlSom <- recommendation[n.feature.interval,n.sample.interval]
    util.info("Recommended SOM size will be used:",env$preferences$dim.1stLvlSom,"x",env$preferences$dim.1stLvlSom)
  }


  # set color schemes
  if (!is.null(env$color.palette.portraits)) # check if given color palette is a valid function
  {
    if( length(environment(env$color.palette.portraits))!=3 || !all( c("colors","ramp") %in% ls(environment(env$color.palette.portraits)) ) )
    {
      util.warn("Invalid value of \"color.palette.portraits\". Using standard scheme")
      env$color.palette.portraits <- colorRampPalette(c("darkblue","blue","lightblue3","green3","yellow2","red2","darkred"))
    }
  } else
  {
    env$color.palette.portraits <- colorRampPalette(c("darkblue","blue","lightblue3","green3","yellow2","red2","darkred"))
  }

  if (!is.null(env$color.palette.heatmaps)) # check if given color palette is a valid function
  {
    if( length(environment(env$color.palette.heatmaps))!=3 || !all( c("colors","ramp") %in% ls(environment(env$color.palette.heatmaps)) ) )
    {
      util.warn("Invalid value of \"color.palette.heatmaps\". Using standard scheme")
      env$color.palette.heatmaps <- colorRampPalette(c("#2C7BB6","#64A4CC","#9CCEE3","#C6E5DB","#ECF6C8","#FEEDAA","#FDC980","#F89D59","#E75B3A","#D7191C"))
    }
  } else
  {
    env$color.palette.heatmaps <- colorRampPalette(c("#2C7BB6","#64A4CC","#9CCEE3","#C6E5DB","#ECF6C8","#FEEDAA","#FDC980","#F89D59","#E75B3A","#D7191C"))
  }

  if(env$preferences$activated.modules$primary.analysis)
  {
    env$files.name <- env$preferences$dataset.name
    while (file.exists(paste(env$files.name, ".RData", sep=""))) {
      env$files.name <- paste(env$files.name, "+", sep="")
    }

    env$output.paths <-
      c("CSV"=paste(env$files.name, "- Results/CSV Sheets"),
        "Summary Sheets Samples"=paste(env$files.name, "- Results/Summary Sheets - Samples") )

    if( !grepl("german",sessionInfo()$locale,ignore.case=TRUE) &&
        !grepl("bioinf.uni-leipzig.de",Sys.info()["nodename"],ignore.case=TRUE) )
      env$csv.function <- write.csv
  }

  env$passedInputChecking <- TRUE
  return(env)
}
