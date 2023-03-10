% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_trackops.R
\name{jadb_addScaledTrack}
\alias{jadb_addScaledTrack}
\title{jadb_addScaledTrack}
\usage{
jadb_addScaledTrack(ids, scale = "zscore", input = "BEADSQ10NU",
  genome = "ce11")
}
\arguments{
\item{IDs}{Vector of JADB ContactExpIDs}
}
\value{
List
}
\description{
jadb_addScaledTrack
}
\examples{
#parallel::mclapply(sprintf('REP\%.3i', 36:42)[-6], addBEADSmapq0TrackZcs), mc.cores = 8)
# sapply(sprintf('REP\%.3i', 3:24), addBEADSmapq0TrackZcs)
}
\seealso{
Other Peaks: \code{\link{addBEADSmapq0TrackZcs}},
  \code{\link{addMapq0TrackZcs}},
  \code{\link{addMapq0Track}},
  \code{\link{addStrandedRNAseq}},
  \code{\link{callPeaksMACS}},
  \code{\link{combinePeaksIDR}},
  \code{\link{combinePeaksToBed}},
  \code{\link{combinePeaks}},
  \code{\link{enrichedRegionsCall}},
  \code{\link{jadb_addTracksFromBAM}},
  \code{\link{meme_chip_local}}
}
\author{
Przemyslaw Stempor
}
