% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNAseq.R
\name{doDiffExpr}
\alias{doDiffExpr}
\title{Run DEseq2 and output results}
\usage{
doDiffExpr(e, design = ~strain)
}
\arguments{
\item{e}{eset}

\item{design}{experiment design}
}
\value{
list
}
\description{
Run DEseq2 and output results
}
\seealso{
Other RNAseq: \code{\link{bam2del}},
  \code{\link{downloadSummarizedEperiment}},
  \code{\link{getAnno}}, \code{\link{getDBdataField}},
  \code{\link{getDEreportForDESeqDataSet}},
  \code{\link{getDEreport}}, \code{\link{getDEseq2res}},
  \code{\link{getStage}}, \code{\link{getStrain}},
  \code{\link{getSummarizedEperiment}},
  \code{\link{ids2SE}}, \code{\link{ids2diffExpr}},
  \code{\link{spikeINsummary}},
  \code{\link{summarizeBAMs}}, \code{\link{wigs2SE}}
}
\author{
Przemyslaw Stempor
}
