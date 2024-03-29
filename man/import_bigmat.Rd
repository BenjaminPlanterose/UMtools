% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_bigmat.R
\name{import_bigmat}
\alias{import_bigmat}
\title{Import big matrices}
\usage{
import_bigmat(filename, nThread = 1)
}
\arguments{
\item{filename}{Name of the file to import}

\item{nThread}{Number of CPU cores to employ}
}
\value{
Matrix
}
\description{
Imports large matrices with the help of data.table::fread
}
\examples{
M = import_bigmat("2021-01-27_M.txt")
}
