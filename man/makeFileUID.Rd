% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/JADBtools.R
\name{makeFileUID}
\alias{makeFileUID}
\title{Make UID}
\usage{
makeFileUID(string = "", prefix = "X")
}
\arguments{
\item{string}{file patch or other string for md5}

\item{prefix}{UID prefix letter, defaults to "X"}
}
\value{
IUD string
}
\description{
Make UID
}
\seealso{
Other dbtools: \code{\link{addGenericFile}},
  \code{\link{correctName}}, \code{\link{formGenericPath}},
  \code{\link{getFilePath}}, \code{\link{getFileUID}}
}
\author{
Przemyslaw Stempor
}
