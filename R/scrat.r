# Creates a new scrat environment
scrat.new <- function(preferences=NULL)
{
  # Init the environment
  env <- new.env()
  env$color.palette.portraits <- NULL
  env$color.palette.heatmaps <- NULL
  env$fdr.g.m <- NULL
  env$files.name <- NULL
  env$gene.info <- NULL
  env$chromosome.list <- NULL
  env$group.silhouette.coef <- NULL
  env$group.colors <- NULL
  env$group.labels <- NULL
  env$gs.def.list <- NULL
  env$samples.GSZ.scores <- NULL
  env$spot.list.correlation <- NULL
  env$spot.list.dmap <- NULL
  env$spot.list.group.overexpression <- NULL
  env$spot.list.kmeans <- NULL
  env$spot.list.overexpression <- NULL
  env$spot.list.samples <- NULL
  env$spot.list.underexpression <- NULL
  env$indata <- NULL
  env$indata.ensID.m <- NULL
  env$indata.gene.mean <- NULL
  env$indata.sample.mean <- NULL
  env$seuratObject <- NULL
  env$metacell.data <- NULL
  env$metadata <- NULL
  env$output.paths <- NULL
  env$pat.labels <- NULL
  env$p.g.m <- NULL
  env$perc.DE.m <- NULL
  env$som.result <- NULL
  env$groupwise.group.colors <- NULL
  env$unique.protein.ids <- NULL
  env$pseudotime.trajectory <- NULL
  env$csv.function <- write.csv2


  # Generate some additional letters
  env$LETTERS <- c(LETTERS, as.vector(sapply(1:10, function(x) {
    return(paste(LETTERS, x, sep=""))
  })))

  env$letters <- c(letters, as.vector(sapply(1:10, function(x) {
    return(paste(letters, x, sep=""))
  })))

  # Set default preferences
  env$preferences <- list(dataset.name = "Unnamed",
													note = "",
                          dim.1stLvlSom = "auto",
                          training.extension = 1,
                          rotate.SOM.portraits = 0,
                          flip.SOM.portraits = FALSE,
                          activated.modules = list( "reporting" = TRUE,
                                                    "primary.analysis" = TRUE,
                                                    "sample.similarity.analysis" = TRUE,
                                                    "geneset.analysis" = TRUE,
                                                    "group.analysis" = TRUE,
                                                    "difference.analysis" = TRUE),
                          database.biomart = "ENSEMBL_MART_ENSEMBL",
                          database.host = "https://jan2020.archive.ensembl.org",
                          database.dataset = "auto",
                          database.id.type = "",
                          standard.spot.modules = "kmeans",
                          spot.coresize.modules = 3,
                          spot.threshold.modules = 0.95,
                          spot.coresize.groupmap = 5,
                          spot.threshold.groupmap = 0.75,
                          adjust.autogroup.number = 0,
                          pseudotime.estimation = TRUE,
                          preprocessing = list(
                            cellcycle.correction = TRUE,
                            create.meta.cells = FALSE,
                            feature.centralization = TRUE,
                            sample.quantile.normalization = TRUE ) )

  # Merge user supplied information
  if (!is.null(preferences))
  {
    env$preferences <-
      modifyList(env$preferences, preferences[names(env$preferences)])
  }
  if(!is.null(preferences$indata))
  {
    env$indata <- preferences$indata
  }
  if(!is.null(preferences$group.labels))
  {
    env$group.labels <- preferences$group.labels
  }
  if(!is.null(preferences$group.colors))
  {
    env$group.colors <- preferences$group.colors
  }

  return(env)
}

# Executes the scrat pipeline.
scrat.run <- function(env)
{
  env$preferences$system.info <- Sys.info()
  env$preferences$session.info <- sessionInfo()
  env$preferences$started <- format(Sys.time(), "%a %d %b %Y %X")

  util.info("Started:", env$preferences$started)
  util.info("Name:", env$preferences$dataset.name)

  #### Preparation & Calculation part ####
  env <- pipeline.checkInputParameters(env)
  if (!env$passedInputChecking) {
    return(env)
  }

  if(env$preferences$activated.modules$reporting)
  {
    # create output directories
    dir.create(paste(env$files.name, "- Results"), showWarnings=FALSE)
    dir.create(paste(env$files.name, "- Results/CSV Sheets"), showWarnings=FALSE)
  }

  if(env$preferences$activated.modules$primary.analysis)
  {
    util.info("Preprocessing Seurat Object")
    env <- pipeline.seuratPreprocessing(env)

    filename <- paste(env$files.name, "pre.RData")
    util.info("Saving environment image:", filename)
    save(env, file=filename)
  }

  if(env$preferences$activated.modules$primary.analysis || env$preferences$activated.modules$geneset.analysis)
  {
    util.info("Loading gene annotation data.")
    env <- pipeline.prepareAnnotation(env)
  }

  if(env$preferences$activated.modules$primary.analysis)
  {
    util.info("Classification of cell cycle phases")
    env <- pipeline.cellcycleProcessing(env)
  }

  if (env$preferences$preprocessing$create.meta.cells)
  {
    env <- pipeline.createMetacells(env)
  }


  if(env$preferences$activated.modules$primary.analysis)
  {
    util.info("Processing SOM. This may take several time until next notification.")
    env <- pipeline.prepareIndata(env)
    env <- pipeline.generateSOM(env)

    filename <- paste(env$files.name, "pre.RData")
    util.info("Saving environment image:", filename)
    save(env, file=filename)

    util.info("Processing Differential Expression Statistics")
    env <- pipeline.calcStatistics(env)

    util.info("Detecting Spots")
    env <- pipeline.detectSpotsSamples(env)
    env <- pipeline.detectSpotsIntegral(env)
    env <- pipeline.patAssignment(env)
    env <- pipeline.groupAssignment(env)
  }

  if (env$preferences$activated.modules$geneset.analysis)
  {
    util.info("Calculating Geneset Enrichment")
    env <- pipeline.genesetStatisticSamples(env)
    env <- pipeline.genesetStatisticIntegral(env)
  }

	 if(!is.null(env$preferences$pseudotime.estimation))
	{
		util.info("Processing Pseudotime Analysis")
		env <- pipeline.pseudotimeEstimation(env)
	}

  if(env$preferences$activated.modules$primary.analysis || env$preferences$activated.modules$geneset.analysis)
  {
    env$indata <- NULL
    env$indata.ensID.m <- NULL

    filename <- paste(env$files.name, ".RData", sep="")
    util.info("Saving environment image:", filename)
    save(env, file=filename)

    if (file.exists(paste(env$files.name, "pre.RData")) && file.exists(filename))
    {
      file.remove(paste(env$files.name, "pre.RData"))
    }
  }

  #### Reporting part ####

  if(env$preferences$activated.modules$reporting)
  {

    pipeline.summarySheetSeurat(env)

    if(ncol(env$metadata) < 1000)
    {
      util.info("Plotting Sample Portraits")
      pipeline.sampleExpressionPortraits(env)
    }

    if ( env$preferences$activated.modules$sample.similarity.analysis && ncol(env$seuratObject) > 2)
    {
      util.info("Plotting (Meta-)cell Similarity Analyses")
      dir.create(file.path(paste(env$files.name, "- Results"), "Sample Similarity Analysis"), showWarnings=FALSE)

      pipeline.sampleSimilarityAnalysisED(env)
      pipeline.sampleSimilarityAnalysisCor(env)
      pipeline.sampleSimilarityAnalysisICA(env)
    }

    util.info("Plotting Summary Sheets (Modules & PATs)")
    pipeline.summarySheetsModules(env)

    if(env$preferences$activated.modules$group.analysis && length(unique(env$group.labels)) >= 2)
    {
      util.info("Processing Group-centered Analyses")
      dir.create(paste(env$files.name, "- Results/Summary Sheets - Groups"), showWarnings=FALSE)
      pipeline.summarySheetsGroups(env)
    }

	  if(!is.null(env$preferences$pseudotime.estimation))
		{
			util.info("Processing Pseudotime Reports")
			pipeline.pseudotimeReport(env)
		}

    util.info("Generating HTML Report")
    pipeline.htmlSummary(env)

  }

  util.info("Finished:", format(Sys.time(), "%a %b %d %X"))

  return(env)
}
