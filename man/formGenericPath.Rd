% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/JADBtools.R
\name{formGenericPath}
\alias{formGenericPath}
\title{Make generic file name}
\usage{
formGenericPath(ContactExpID, EXTABLE = "labrnaseq", Processing = NA,
  Resolution = NA, Scale = NA)
}
\arguments{
\item{ContactExpID}{experiment ID string}

\item{EXTABLE}{experiment table in db schema}

\item{scale}{file name parameter}

\item{processing}{file name parameter}

\item{format}{file name parameter}
}
\value{
name string
}
\description{
Make generic file name
}
\seealso{
Other dbtools: \code{\link{addGenericFile}},
  \code{\link{correctName}}, \code{\link{getFilePath}},
  \code{\link{getFileUID}}, \code{\link{makeFileUID}}
}
\author{
Przemyslaw Stempor
}
