\name{generate.dag.data}
\alias{generate.dag.data}
\alias{rdag}
\alias{moralize}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generate nonlinear data from DAGs
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
These funtions create distributions on directed acyclic graphs. \code{rdag} generates a random DAG with a given number of edges by selecting 'nedges' at random from (p choose 2) possible edges, 
\code{moralize} gets the moral graph from a DAG, and
\code{generate.dag.data} generates non-linear data by assigning each edge a cubic polynomial basis with random coefficients.
}
\usage{
rdag(p, nedges)
moralize(g)
generate.dag.data(g, n, basesd = 1, basemean = 0, bfuns = function(x){cbind(x, x^2, x^3)}, 
	funclist = NULL, usenorm = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
%%     ~~Describe \code{g} here~~
a directed graph, as an 'igraph' object
}
  \item{p}{
%%     ~~Describe \code{n} here~~
number of vertices
}
  \item{n}{
%%     ~~Describe \code{n} here~~
number of observations
}
  \item{nedges}{
%%     ~~Describe \code{n} here~~
number of edges
}
\item{bfuns}{
%%     ~~Describe \code{basemean} here~~
the basis functions for the structural equations. Note that when the basis functions generated, they are centered and scaled to have variance 1, so similar coefficients correspond to similar amounts of variance explained. Ignored if `funclist` is supplied.
}
  \item{basesd}{
%%     ~~Describe \code{basesd} here~~
standard deviation of the random coefficients assigned to the basis functions. Ignored if `funclist` is supplied.
}
  \item{basemean}{
%%     ~~Describe \code{basemean} here~~
means of the (random) coefficients assigned to the basis functions. Ignored if `funclist` is supplied.
}
\item{funclist}{
%%     ~~Describe \code{basemean} here~~
p by p list of functions determining the structural equations. funclist[[i]][[j]](x) is the effect of feature j on feature i. If this is omitted, the functions are generated at random from the bases supplied in bfuns.
}
  \item{usenorm}{
%%     ~~Describe \code{usenorm} here~~
logical. whether to use normal or uniform errors
}
}
\details{
Multivariate distributions with complicated conditional dependence structures corresponding to a particular graph are difficult to construct in general. However, constructing complicated distributions from a DAG is straighforward. These functions are meant to facilitate construction of complicated distributions on a DAG, and obtain the coresponding conditional independence structure. 

generate.dag.data generates Normal(basemean,basesd) coefficients for the basis functions coresponding to each edge. This gives a function f_ij(x_j), which is standardized to have variance 1. Then data is generated for a feature conditioned on its parents by

x_i = sum_(j in parents(i)) f_ij + noise 

where `noise' is by default N(0,1). These data are returned, along with the generated functions which can be used in subsequent simulations.
}
\value{
`rdag' returns an `igraph` object. 

`generate.dat.data' returns a list with two elements: an n by p matrix 'X' and a p by p list `funclist', where funclist[[i]][[j]] is the effect of feature j on feature i.

`moralize` returns an undirected igraph object.
}
\author{
%%  ~~who you are~~
Arend Voorman
}
\references{
Voorman, Shojaie and Witten (2013). Graph Estimation with Joint Additive Models. Submitted to Biometrika. available on ArXiv or from authors upon request
}
\seealso{
\code{\link{igraph}}
\code{\link{spacejam}}
\code{\link{SJ}}
\code{\link{plot.SJ}}
}
\examples{
#########Create graph and distribution used in Figure 2 of Voorman, Shojaie and Witten (2013):
p <- 100 #variables
n <- 50 #observations

#Generate Graph
set.seed(20)
g <- rdag(p,80)
mylayout <- layout.fruchterman.reingold(g)

par(mfrow=c(1,2))
plot(g, layout = mylayout, edge.color = "gray50", 
        vertex.color = "red", vertex.size = 3, vertex.label = NA, 
        edge.arrow.size = 0.4)
plot(moralize(g), layout = mylayout, edge.color = "gray50", 
        vertex.color = "red", vertex.size = 3, vertex.label = NA, 
        edge.arrow.size = 0.4)

#create a distribution on the DAG using cubic polynomials with random normal coefficients 
#with standard deviations of 1, 0.5 and 0.5, (i.e. giving more weight to linear association than quadratic or cubic)
data <- generate.dag.data(g,n,basesd=c(1,0.5,0.5))
X <- data$X

#Fit conditional independence graph at one lambda 
fit1 <- SJ(X, lambda = 0.6)

###For additional replications using the same DAG distribution use e.g.
data <- generate.dag.data(g,n,funclist = data$funclist)

}