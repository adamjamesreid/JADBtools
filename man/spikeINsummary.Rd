% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNAseq.R
\name{spikeINsummary}
\alias{spikeINsummary}
\title{Get spikeIN summary from IDs}
\usage{
spikeINsummary(ContactExpIDs = paste0("rFB0", 26:32))
}
\arguments{
\item{ContactExpID}{vector of IDs}
}
\value{
SummarizedEperiment
}
\description{
Get spikeIN summary from IDs
}
\seealso{
Other RNAseq: \code{\link{bam2del}},
  \code{\link{doDiffExpr}},
  \code{\link{downloadSummarizedEperiment}},
  \code{\link{getAnno}}, \code{\link{getDBdataField}},
  \code{\link{getDEreportForDESeqDataSet}},
  \code{\link{getDEreport}}, \code{\link{getDEseq2res}},
  \code{\link{getStage}}, \code{\link{getStrain}},
  \code{\link{getSummarizedEperiment}},
  \code{\link{ids2SE}}, \code{\link{ids2diffExpr}},
  \code{\link{summarizeBAMs}}, \code{\link{wigs2SE}}
}
\author{
Przemyslaw Stempor
}
