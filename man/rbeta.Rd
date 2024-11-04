\name{rbeta}
\alias{rbeta}
\title{
 The Beta Random Number Generating Function
}
\description{
 Random generation for the beta distribution with parameters \code{shape1} and
 \code{shape2}.
}
\usage{
 rbeta(n, shape1, shape2)
}
\arguments{
 \item{n}{Number of beta random numbers to generate. If \code{length(n) > 1}, the length 
    is taken to be the number required.}
 \item{shape1, shape2}{Positive shape parameters.}
}
\details{
 The beta distribution with parameters \code{shape1}\eqn{ = a} and \code{shape2}\eqn{ = b} 
 has density
 \deqn{\frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)} x^{a-1} (1-x)^{b-1}}{%
 \Gamma(a+b)/(\Gamma(a)\Gamma(b)) x^(a-1)(1-x)^(b-1)}
 for \eqn{a > 0, b > 0} and \eqn{0 \le x \le 1}.
 
 The mean is \eqn{\frac{a}{a+b}}{a/(a+b)} and the variance is
 \eqn{\frac{ab}{(a+b)^2 (a+b+1)}}{ab/((a+b)^2 (a+b+1))}.
 
 \code{rbeta} basically utilizes the following guideline primarily proposed by Hung 
 \emph{et al.} (2009) for generating beta random numbers.
 \itemize{
  \item When \eqn{max(}\code{shape1}\eqn{, }\code{shape2}\eqn{) < 1}, the B00 algorithm 
     (Sakasegawa, 1983) is used;
  \item When \code{shape1}\eqn{ < 1 < }\code{shape2} or \code{shape1}\eqn{ > 1 > }\code{shape2}, 
     the B01 algorithm (Sakasegawa, 1983) is used;
  \item When \eqn{min(}\code{shape1}\eqn{, }\code{shape1}\eqn{) > 1}, the B4PE algorithm
     (Schmeiser and Babu, 1980) is used if one papameter is close to 1 and the other is large
     (say \eqn{> 4}); otherwise, the BPRS algorithm (Zechner and Stadlober, 1993) is used.
 }
}
\value{
 \code{rbeta} generates beta random numbers.
}
\source{
 \code{rbeta} uses a C translation of
 
 Y. C. Hung and N. Balakrishnan and Y. T. Lin (2009),
  Evaluation of beta generation algorithms,
  \emph{Communications in Statistics - Simulation and Computation},
  \bold{38}:750--770.
}
\references{
 Y. C. Hung and N. Balakrishnan and Y. T. Lin (2009),
  Evaluation of beta generation algorithms,
  \emph{Communications in Statistics - Simulation and Computation}, \bold{38}, 750--770.
  
 H. Sakasegawa (1983),
  Stratified rejection and squeeze method for generating beta random numbers,
  \emph{Annals of the Institute Statistical Mathematics}, \bold{35}, 291--302.
  
 B.W. Schmeiser and A.J.G. Babu (1980),
  Beta variate generation via exponential majorizing functions,
  \emph{Operations Research}, \bold{28}, 917--926.
  
 H. Zechner and E. Stadlober (1993),
  Generating beta variates via patchwork rejection,
  \emph{Computing}, \bold{50}, 1--18.
}
\author{
 Ching-Wei Cheng <aks43725@gmail.com>,\cr
 Ying-Chao Hung <hungy@nccu.edu.tw>,\cr
 Narayanaswamy Balakrishnan <bala@univmail.cis.mcmaster.ca>
}
\seealso{
 \code{\link[stats]{rbeta}} in package \pkg{stats}.
}
\examples{
 library(rBeta2009)
 rbeta(10, 0.7, 1.5)
}
\keyword{beta}
