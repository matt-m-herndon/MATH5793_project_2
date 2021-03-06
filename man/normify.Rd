% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normal.R
\name{normify}
\alias{normify}
\title{Run normal on a dataset and return a "normified" dataset, transformed when necessary}
\usage{
normify(data, ...)
}
\arguments{
\item{data}{2D Dataframe to normify.}
}
\description{
Run normal on a dataset and return a "normified" dataset, transformed when necessary
}
\examples{

# Load iris dataset 
data('iris')

normifed_data = normify(iris[1:25,], significance=5)

}
