\name{scrat.genesets}
\alias{scrat.genesets}
\title{Additional literature genesets}
\description{Genesets collected from publications and independent analyses.}
\usage{
data(scrat.genesets)
}
\format{
The data set is stored in RData (binary) format. Each element of the list represents one distinct gene set and contains the Ensembl-IDs of the member genes.
}
\details{
The scrat package allows for analysing the biological background of the samples using predefined sets of genes of known biological context. A large and diverse collection of such gene sets is automatically derived from the Gene Ontology (GO) annotation database using biomaRt interface. \code{scrat.genesets} contains about 150 additional gene sets collected from literature on cellular systems, and from previous analyses using the SOM pipeline.
}