\name{scrat.new}
\alias{scrat.new}
\title{
  Initialize the scrat pipeline.
}
\description{
  This function initializes the scrat environment and sets the preferences.
}
\usage{
  scrat.new(preferences)
}
\arguments{
  \item{preferences}{
    list with the following parameters:

    \itemize{
      \item{indata}: input data matrix containing the expression values or an \code{scater::SCESet} object (see 'Details' and 'Examples')
      \item{group.labels}: cell assignment to a distinct group, celltype or class (character; "auto" or one label for each cell; may be given with \code{indata SCESet})
      \item{group.colors}: colors of the cells for visualizations (character; one color for each cell; may be given with \code{indata SCESet})
      \item{dim.1stLvlSom}: dimension of the SOM grid (integer, >5); use "auto" to use size recommendation
      
      \item{database.dataset}: type of ensemble dataset queried using biomaRt interface (character); use "auto" to detect parameter automatically; for details see vignette 
      \item{database.id.type}: type of rowname identifier in biomaRt database (character); obsolete if \code{database.dataset="auto"} 

      \item{standard.spot.modules}: spot modules utilized in downstream analyses such as PAT detection and unsupervised cell grouping (character, one of \{"overexpression", "kmeans", "group.overexpression", "dmap"\})
    
      \item {preprocessing}: activates/deactivates preprocessing steps (list):
        \itemize{
          \item{count.processing}: enables count filtering and normalization when input data are read counts (boolean).
          \item{cellcycle.correction}: enables removal of cell cycle effect from expression values (only available for human and mouse organisms; boolean).    
        	\item{feature.centralization}: enables or disables
            centralization of the features (boolean).
        	\item{sample.quantile.normalization}: enables quantile
            normalization of the cells (boolean).
        }

      \item{pseudotime.estimation}: parameters for wanderlust pseudotime estimation: NULL deactivates pseudotime analysis, TRUE uses standard parameter setting, or a list:
        \itemize{
      		\item{n.waypoints}: number of waypoint cells (integer, >2 )
      		\item{n.iterations}: number of iterations (integer, >1 )
      		\item{k}: number of neighbors of each node in graph (integer, >2 )
      		\item{I}: number of neighbors cells (integer,  1<I<k )
      		\item{initiator.sample}: index or name of initiator cell (integer or character)
        }

    }
  }
}
\details{
scrat requires input of mapped read counts or preprocessed expression data. Additionally, it is recommended to provide group information about the cells (e.g. celltypes), otherwise groups will be assigned automatically based on cells' expression landscapes.
\code{indata} is to be provided as \code{scater::SCESet} object, or as two-dimensional numerical matrix (columns and rows represent the cells and genes, respectively). scrat can apply read count preprocessing or use preprocessed expression values.
Please check the vignette for more details on the parameters.
}
\value{
  A new scrat environment which is passed to \code{scrat.run}.
}
\examples{
env <- scrat.new(list(dataset.name="Example",
                        dim.1stLvlSom="auto",
                        database.dataset="auto",
                        standard.spot.modules="kmeans",
                        
                        preprocessing = list(
                          count.processing = FALSE,
                          cellcycle.correction = FALSE,
                          feature.centralization = TRUE,
                          sample.quantile.normalization = TRUE ), 
                        
                        pseudotime.estimation = list( 
                          n.waypoints=20,
                          n.iterations=20,
                          k=30,
                          I=5,
                          initiator.sample=1) ) )


# definition of indata, group.labels and group.colors
env$indata = matrix( runif(1000), 100, 10 )
env$group.labels = c( rep("class 1", 5), rep("class 2", 4), "class 3" )
env$group.colors = c( rep("red", 5), rep("blue", 4), "green" )

}
