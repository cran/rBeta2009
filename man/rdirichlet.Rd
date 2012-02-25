\name{rdirichlet}
\alias{rdirichlet}
\title{
 The Dirichlet Random Vector Generating Function
}
\description{
 The function to generate random vectors from the Dirichlet distribution.
}
\usage{
 rdirichlet(n, shape)
}
\arguments{
 \item{n}{Number of Dirichlet random vectors to generate. If \code{length(n) > 1},
	the length is taken to be the number required.}
 \item{shape}{Vector with \code{length(shape) >= 2} containing positive shape parameters of 
  the Dirichlet distribution. If \code{length(shape) = 2}, it reduces to the beta 
  generating function.}
}
\details{
 The Dirichlet distribution is the multidimensional generalization of the beta 
 distribution.
 
 A \eqn{k}-variate Dirichlet random vector \eqn{(x_1,\ldots,x_k)}{(x[1],\ldots,x[k])} has 
 the joint probability density function
 \deqn{\frac{\Gamma(\alpha_1+\dots+\alpha_{k+1})}{\Gamma(\alpha_1)\dots\Gamma(\alpha_{k+1})}
 x_1^{\alpha_1-1}\dots x_k^{\alpha_k-1}\left(1-\sum_{i=1}^k x_i\right)^{\alpha_{k+1}-1},}{%
 \Gamma(\alpha[1]+\dots+\alpha[k+1])(\Gamma(\alpha[1])\dots\Gamma(\alpha[k+1]))
 x[1]^(\alpha[1]-1)\dots x_k^(\alpha[k]-1)(1-\sum_{i=1}^k x[i])^(\alpha[k+1]-1),}
 where \eqn{x_i \ge 0}{x[i] \ge 0} for all \eqn{i = 1, \ldots, k},
 \eqn{\sum_{i=1}^k x_i \leq 1}{\sum_{i=1}^k x[i] \le 1}, and 
 \eqn{\alpha_1, \ldots, \alpha_{k+1}}{\alpha[1], \ldots, \alpha[k+1]} are positive shape
 parameters.
 
 \code{rdirichlet} generates the Dirichlet random vector by utilizing the transformation
 method based on beta variates and three guidelines introduced by Hung \emph{et al.} (2011).
 The three guidelines include: how to choose the fastest beta generation algorithm, how to
 best re-order the shape parameters, and how to reduce the amount of arithmetic operations.
}
\value{
 \code{rdirichlet()} returns a matrix with \code{n} rows, each containing a single Dirichlet
 random vector.
}
\source{
 \code{rdirichlet} uses a C translation of
 
 Y. C. Hung and N. Balakrishnan and C. W. Cheng (2011),
  Evaluation of algorithms for generating Dirichlet random vectors,
  \emph{Journal of Statistical Computation and Simulation}, \bold{81}, 445--459.
}
\references{
 Y. C. Hung and N. Balakrishnan and C. W. Cheng (2011),
  Evaluation of algorithms for generating Dirichlet random vectors,
  \emph{Journal of Statistical Computation and Simulation}, \bold{81}, 445--459.
}
\author{
 Ching-Wei Cheng <aks43725@gmail.com>,\cr
 Ying-Chao Hung <hungy@nccu.edu.tw>,\cr
 Narayanaswamy Balakrishnan <bala@univmail.cis.mcmaster.ca>
}
\seealso{
 \code{\link[MCMCpack]{rdirichlet}} in package \pkg{MCMCpack}.\cr
 \code{\link[gtools]{rdirichlet}} in package \pkg{gtools}.
}
\examples{
library(rBeta2009)
rdirichlet(10, c(1.5, 0.7, 5.2, 3.4))

}
\keyword{Dirichlet}
