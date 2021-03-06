% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{vect2Mat}
\alias{vect2Mat}
\title{vect2Mat : transform a vector into a  matrix  (replicating symetric elements if needed)}
\usage{
vect2Mat(V, symmetric, diag = FALSE)
}
\arguments{
\item{V}{: vector}

\item{symmetric}{: TRUE is the network is not  directed. FALSE otherwise}

\item{diag}{: FALSE if the self interactions are neglected. TRUE otherwise (default value =  FALSE)}
}
\value{
matrix fulfilled with the elements of the matrix
}
\description{
vect2Mat : transform a vector into a  matrix  (replicating symetric elements if needed)
}
\examples{
V <- c(0,15);
M <- vect2Mat(V, symmetric = TRUE, diag=FALSE)
}
