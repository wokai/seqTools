\name{kMerIndex}
%
\alias{kMerIndex}
%
\title{
kMerIndex function: Returns array index for given DNA k-mers.
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{
For each k, there exist $4^k$ DNA k-mers. Many functions inside this package
    return values where DNA k-mers appear as array indices. \code{kMerIndex} 
    can be used for extraction of count values for special k-mers by provision
    of index values.
}
%
\usage{
kMerIndex(kMers, k=nchar(kMers)[1], base=1)
}
%
\arguments{
    \item{kMers}{\code{character}. Vector of equal sized character strings.
        The number of characters in each string must be =k 
        (i.e. \code{all(nchar(kMers)==k)})}
    \item{k}{\code{integer}. Length of k-mer.}
    \item{base}{\code{integer}. Value must be 0 or 1 (i.e. length(base)==1). 
        For \code{base=0} the returned index is 0-based (i.e. the index of the
        first k-mer (AAA..)) is 0. Otherwise the index is 1-based.}
}
%
\details{Maximal allowed value for k is 12.}
%
\value{\code{integer}.}
%
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
kMerIndex(c("AACC", "ATAA"))
kMerIndex(c("AA","AC"), base=1)
kMerIndex(c("AA","AC"), base=0)
}
\keyword{kMerIndex}
