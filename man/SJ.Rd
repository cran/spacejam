\name{SJ}
\alias{SJ}
\alias{SJ.dag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Estimate a graph with Spacejam
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
These functions estimate graphs using flexible node-wise regressions, employing a group-lasso penalty to encourage sparsity. Details are given in Voorman, Shojaie and Witten (2013).
}
\usage{
SJ(X, bfun = bs, lambda=NULL, length =NULL, verbose = FALSE, b0 = NULL, maxit = 100, 
	tol = .Machine$double.eps^0.25,G.max = NULL)
	
SJ.dag(X, bfun = bs, lambda=NULL, length =NULL, ord = 1:p, verbose = FALSE, b0 = NULL, 
	maxit = 100, tol = .Machine$double.eps^0.25)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
an n by p matrix of observations
}
  \item{bfun}{
%%     ~~Describe \code{bfun} here~~
a function to generate bases for each of the p features (e.g. 'bs' or \code{function(x){cbind(x,x^2)}}). The default is \code{bs}, which generates cubic polynomials (i.e. spacejam in 3 dimensions)
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
penalty terms. These should be between 0 and k, where k is the number of basis functions. If none are specified, the `\code{length}' argument specifies the number of penalty terms to choose automatically (default is 10). 
}
  \item{length}{
%%     ~~Describe \code{lambda} here~~
number of penalty terms. This is ignored if lambda is specified. 
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
logical.  whether to print progress indicator.
}
  \item{b0}{
%%     ~~Describe \code{b0} here~~
a p*k by p matrix of coefficients to be used as initial values. Typically obtained from the `betas` value of an SJ object.
}
  \item{maxit}{
%%     ~~Describe \code{maxit} here~~
maximum number of iterations
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
Tolerance for convergence.  This is the maximum change in coefficients in subsequent iterations.
}
\item{ord}{
%%     ~~Describe \code{tol} here~~
for SJ.dag: The causal ordering of the variables, according to the columns of X. If none is given, assumed to be 1:ncol(X) i.e. the first column of X is first in the causal ordering, and the final column is last in the causal ordering.
}
\item{G.max}{
	An adjacency matrix of dimension ncol(X) by ncol(X) specifying all potential edges to be considered. This is a Useful for screening rules, where many edges are not considered, or incorporating known structures of the graph.
}
}
\details{
This implements the method described in Voorman, Shojaie and Witten `Graph estimation with joint additive models`, which regresses each feature on the others using generalized additive models, subject to sparsity inducing penalties.

The function is designed to estimate the graph for a sequence of tuning parameters, using warm starts to improve speed. In addition, an `active set' approach is used as described in Friedman et al (2010). 
}
\value{
an object of class 'SJ' or 'SJ.dag'.

Among some internal variables, these objects include the elements 

\item{graphs}{
	 a list of the estimates, each of class igraph
}
\item{lambda}{
	 a vector of the corresponding tuning parameters
}
\item{bic}{ 
	a vector of estimated BIC criteria.
}
\item{dfs}{
	 a p x length matrix of the estimated degrees of freedom for each feature and tuning parameter.
}
\item{rss}{
	a p x length matrix of the residual sum of squares for each feature and tuning parameter.
}
}
\references{
Voorman, Shojaie and Witten (2013). Graph Estimation with Joint Additive Models. Submitted to Biometrika. available on ArXiv or from authors upon request

Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of
  Statistical Software, 33(1), 1-22. URL http://www.jstatsoft.org/v33/i01/.
}
\author{
%%  ~~who you are~~
Arend Voorman
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{spacejam}}
\code{\link{generate.dag.data}}
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

#Fit conditional independence graph at one lambda , using the default basis functions (cubic polynomials).
fit1 <- SJ(X, lambda = 0.6)

#Fit conditional independence graph at 10 (hopefully reasonable) lambdas:
fit2 <- SJ(X, length = 10)

#Fit conditional independence graph using quadratic basis functions:
fit3 <- SJ(X, bfun = function(x){cbind(x,x^2)}, length = 10)

#Fit the DAG using default causal ordering 1:p, and at 10 lambdas
fit4 <- SJ.dag(X, length = 10)
fit4

#plot the DAGs, and the true graph
par(mfrow=c(1,3))
plot(g, layout=mylayout, edge.color = "gray50", vertex.color = "red", vertex.size = 3, vertex.label = NA, edge.arrow.size = 0.4, main="True DAG")
plot(fit4, layout = mylayout, which= 4, main= paste0("lambda = ",round(fit4$lambda[4],2) ))
plot(fit4, layout = mylayout, main = "min BIC")

###For additional replications using the same DAG distribution use e.g.
data <- generate.dag.data(g,n,funclist = data$funclist)


#### Screen out edges whose corresponding nodes have low spearman correlation:
# useful for approximating spacejam in high dimesnions
S <- cor(data$X, method = "spearman")
G.max <- S > 0.1
mean(G.max) #~75\% of edges removed

system.time({fit.screen <- SJ(X, G.max = G.max,length=20)})
system.time({fit.noscreen <- SJ(X,length = 20)})

plot(fit.screen, layout=mylayout)
plot(fit.noscreen, layout=mylayout)
}
\keyword{ graphs }
\keyword{ optimize}