% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filesStats.R
\name{bamStats}
\alias{bamStats}
\title{Returms bam file stats}
\usage{
bamStats(f, aln = NULL)
}
\arguments{
\item{f}{URL or path to BAM file}
}
\value{
stats string
}
\description{
String is formated as: 
all=%d, aligned=%d[%.2f%%], mapq10=%d[%.2f%%], unique10=%d[%.2f%%]
}
\author{
Przemyslaw Stempor
}
