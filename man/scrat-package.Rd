\name{scrat-package}
\alias{scrat-package}
\alias{scrat}
\docType{package}
\title{
Single Cell - R Analysis Toolbox
}
\description{
This package translates single cell expression data into metadata of reduced
dimension. It provides various sample-centered and group-centered visualizations,
sample similarity analyses and pseudotime analysis. The underlying
SOM algorithm combines feature clustering, multidimensional scaling and dimension
reduction, along with strong visualization capabilities. It enables extraction
and description of functional expression modules inherent in the data.
The results are given within a separate folder and can be browsed using the
summary HTML file.
}
\details{
\tabular{ll}{
Package: \tab scrat\cr
Type: \tab Package\cr
Version: \tab 2.0.2\cr
Date: \tab 2020-01-17\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Author: Henry Loeffler-Wirth <wirth@izbi.uni-leipzig.de> and Gray Camp <gray_camp@eva.mpg.de>

Maintainer: Henry Loeffler-Wirth <wirth@izbi.uni-leipzig.de>
}
\references{
Loeffler-Wirth, Kalcher, Binder: oposSOM: R-package for high-dimensional portraying %
of genome-wide expression landscapes on Bioconductor. (Bioinformatics 2015)

Camp, G., Loeffler-Wirth, H. et al.: Multilineage communication regulates human liver %
bud development from pluripotency. (Nature 2017)
}
\keyword{ package }
\examples{
## Example with artificial data
env <- scrat.new(list(dataset.name="Example",
                        dim.1stLvlSom=10))

env$indata <- matrix(runif(10000), 1000, 10)

scrat.run(env)

## Real data example - This will take several minutes
# env <- scrat.new(list(dataset.name="Macrophages Showcase"))
#
# data(macro.SCESet.counts)
# env$indata <- macro.SCESet.counts
#
# scrat.run(env)
}
