\name{ascii2char}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Alias
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\alias{ascii2char}
\alias{char2ascii}
%
\title{ascii2char: Converting ASCII encoded values to character values.}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Description
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\description{
\code{ascii2char} calculates character representations for given phred values.
    \code{char2ascii} returns phred values for given ASCII encoded 
    representations (the reverse transformation of \code{ascii2char}).}
\usage{
ascii2char(x, multiple=FALSE)
char2ascii(c)
}
\arguments{
    \item{x}{\code{numeric}. Vector with ASCII values. All values must be in
        1:255. Other values produce an error.}
    \item{multiple}{\code{logical}. For \code{'FALSE'} (the default), all 
        characters are combined into one single string (i.e. a character vector
        of length 1). For \code{'TRUE'}, single characters are combined into
        a vector.}
    \item{c}{\code{character}. Vector of length 1 (Longer vectors will 
        generate Warnings).}
}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Details
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\details{The functions are only wrappers for convenience. \code{char2ascii} 
    is defined as \code{strtoi(charToRaw(c), base = 16L)}. \code{ascii2char} 
    is defined as \code{rawToChar(as.raw(x), multiple)}.}
\value{\code{ascii2char} returns \code{character}. \code{char2ascii} 
    returns \code{integer}.}
\references{
Ewing B, Green P
Base-Calling of Automated Sequencer Traces Using Phred. II. Error Probabilities
Genome Research 1998 8(3): 186-194
}
\author{Wolfgang Kaisers}
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% Examples
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
\examples{
ascii2char(97:101, multiple=FALSE)
ascii2char(97:101, multiple=TRUE)
char2ascii("abcde")
char2ascii(paste("a", "b", "c", collapse=""))
ascii2char(char2ascii("abcde"))
}
\seealso{getPhredTable}
\keyword{ascii2char}
