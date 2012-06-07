\name{plot.SJ}
\alias{plot.SJ}
\alias{plot.SJ.dag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
plot an object of class \code{SJ} or \code{SJ.dag}
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function plots an SJ or SJ.dag 
}
\usage{
\method{plot}{SJ}(x, which = NULL, layout = NULL, \dots) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{X} here~~
an object of class SJ or SJ.dag
}
  \item{which}{
%%     ~~Describe \code{bfun} here~~
The index of which to plot. If not specified, the penalty with the smallest BIC is used.
}
  \item{layout}{
%%     ~~Describe \code{lambda} here~~
the layout of the graph to use.  If not specified, \code{\link{layout.fruchterman.reingold}} is used. 
}
  \item{\dots}{
additional parameters to be passed to \code{\link{plot.igraph}}.
}
}
\details{
This function plots a graph in SJ or SJ.dags 'graph' field, using some sane degaults for vertex size and color.
}
\value{
returns the layout used, invisibly.
}

\references{
Voorman, Shojaie and Witten (2013). Graph Estimation with Joint Additive Models. Submitted to Biometrika. available on ArXiv or from authors upon request
}
\author{
%%  ~~who you are~~
Arend Voorman
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{plot.igraph}}
\code{\link{spacejam}}
\code{\link{SJ}}
\code{\link{generate.dag.data}}
}
\examples{
p <- 100 #variables
n <- 50 #observations

#Generate Data
set.seed(20)
g <- rdag(p,80)
data <- generate.dag.data(g,n,basesd=c(1,0.5,0.5))
X <- data$X

#Fit conditional independence graph for sequence of 10 lambdas
fit1 <- SJ(X, length = 10)

par(mfrow=c(1,2))
layout <- plot(fit1, main = "min BIC")
plot(fit1, which=5, layout = layout, main = paste0("lambda = ",round(fit1$lambda[5],3)))
}