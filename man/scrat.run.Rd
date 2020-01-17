\name{scrat.run}
\alias{scrat.run}
\title{
Execute the scrat pipeline.
}
\description{
This function realizes the complete pipeline functionality: single gene
expression values are culstered to metagenes using a self-organizing map. Based
on these metagenes, visualizations (e.g. expression portraits), downstreaming
sample similarity analyses (e.g. hierarchical clustering, ICA) and pseudotime analysis are performed.
The results are given within a separate folder and can be browsed using the
summary HTML file.
}
\usage{
scrat.run(env)
}
\arguments{
  \item{env}{
    the scrat environment created with \code{scrat.new} according to the users' preferences
  }
}

\examples{
# Example with artificial data
env <- scrat.new(list(dataset.name="Example",
                        dim.1stLvlSom=10))

env$indata <- matrix(runif(1000), 100, 10)

scrat.run(env)

}
