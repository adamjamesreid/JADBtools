% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/JADBtools.R
\name{getFileUID}
\alias{getFileUID}
\title{Get file UID}
\usage{
getFileUID(ID, format = ".", processing = ".", scale = ".", eq = FALSE)
}
\arguments{
\item{ID}{ContactExperimetID}

\item{format}{filter regex}

\item{processing}{filter regex}

\item{scale}{filter regex}

\item{eq}{t}

\item{url}{return url or just the path to the file}
}
\value{
URL string
}
\description{
Accepts regular expressions.
}
\seealso{
Other dbtools: \code{\link{addGenericFile}},
  \code{\link{correctName}}, \code{\link{formGenericPath}},
  \code{\link{getFilePath}}, \code{\link{makeFileUID}}
}
\author{
Przemyslaw Stempor
}
