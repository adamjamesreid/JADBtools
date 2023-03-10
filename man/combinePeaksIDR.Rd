% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peakTools.R
\name{combinePeaksIDR}
\alias{combinePeaksIDR}
\title{Combines peak in replicates and exports to bed file}
\usage{
combinePeaksIDR(ids, mode = "union")
}
\arguments{
\item{mode}{union or intersection}

\item{IDs}{Vector of JADB ContactExpIDs}
}
\value{
bed file path
}
\description{
Combines peak in replicates and exports to bed file
}
\examples{
#combinePeaksToBed(IDs)
}
\seealso{
Other Peaks: \code{\link{addBEADSmapq0TrackZcs}},
  \code{\link{addMapq0TrackZcs}},
  \code{\link{addMapq0Track}},
  \code{\link{addStrandedRNAseq}},
  \code{\link{callPeaksMACS}},
  \code{\link{combinePeaksToBed}},
  \code{\link{combinePeaks}},
  \code{\link{enrichedRegionsCall}},
  \code{\link{jadb_addScaledTrack}},
  \code{\link{jadb_addTracksFromBAM}},
  \code{\link{meme_chip_local}}
}
\author{
Przemyslaw Stempor
}
