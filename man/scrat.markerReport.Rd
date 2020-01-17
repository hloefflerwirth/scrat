\name{scrat.markerReport}
\alias{scrat.markerReport}
\title{
Generate a report PDF for a given set of marker genes.
}
\description{
This function generates expression profile plots and mapping into SOM coordinate system for a given set of marker genes.
The resulting file can be found in "Marker Set Reports" within the results folder.
}
\usage{
scrat.markerReport(env,filepath)
}
\arguments{
  \item{env}{
    the scrat environment after \code{scrat.run}.
  }
  \item{filepath}{
    the path of a txt file containing the marker genes (one gene per line).
  }
}