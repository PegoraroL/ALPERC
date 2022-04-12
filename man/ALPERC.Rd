\name{ALPERC}
\alias{ALPERC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ALPERC
}
\description{
ALPERC is a model-agnostic Active Learning (AL) algorithm for noisy physical experiments with multiple responses based on nonparametric ranking and clustering.
}
\usage{
ALPERC(
n_add,
D_cand,
sigma_cand,
S,
strategy,
varimp_distance,
n_clust,
n_boot,
alpha_rank,
seed_rank,
paral
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n_add}{The size of the batch to add at the AL iteration. An integer.
}
  \item{D_cand}{A data frame including the candidate experimental configurations for AL. The first column is a factor column that contains the identifier (ID) of the candidates, the following columns include the factor levels. The column headers (except from the first) are the names of the factors.
}
  \item{sigma_cand}{A data frame including the predictive uncertainty of the candidate points in D_cand. The first column is a factor column that contains the identifier (ID) of the candidates, the following columns include the uncertainty of predictions of the candidate configurations for each response (>2). The column headers (except from the first) are the names of the responses.
  }
  \item{S}{A data frame including the variable importance of each factor. The first column is a character column that contains the names of the factors (equal to the column headers in D_cand), the following columns include the variable importance of the factors corresponding to each response (>2). The column headers (except from the first) are the names of the responses. If variable importance is not available, set S=NULL.
}
  \item{strategy}{It is either "exploration" (only unique candidates) or "exploitation" (can include replicates).
}
  \item{varimp_distance}{It is either TRUE or FALSE. In the first case, the assignment of each candidate configuration to a cluster is adjusted by the importance of each dimension. If varimp_distance=FALSE, S is not used.
}
  \item{n_clust}{The number of clusters, it can either be NULL or an integer. If n_clust=NULL, the best number of clusters is chosen with the Silhouette method, otherwise it is set equal to the number indicated.
}
  \item{n_boot}{The number of bootstrap iterations for the permutations used for the nonparametric ranking procedure.
}
  \item{alpha_rank}{The significance level for the computation of the nonparametric ranking procedure.
}
  \item{seed_rank}{A seed that is used to initialize the permutations.
}
  \item{paral}{Either TRUE or FALSE. It indicates whether parallel computation should be employed (default is FALSE). If paral=TRUE, a cluster must be initialized (e.g. via "doSNOW" library).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{strategy}{Is the strategy chosen for ALPERC implementation.}
\item{seed_rank}{Is the seed chosen to initialize the nonparametric ranking procedure.}
\item{varimp_distance}{Indicates whether the variable importance is used or not.}
\item{clust_choice}{A tibble that includes the number of clusters and corresponding Silhouette index.}
\item{nclust_criterion}{Indicates the strategy used for the selection of the number of clusters.}
\item{best_nclust}{The number of clusters considered.}
\item{LAMBDA}{A tibble that includes the ID of each candidate, the factor levels combination, the position in the rank and the cluster to which each configuration is assigned.}
\item{D_add_j}{A tibble that includes the n_add configurations from D_add that should be added according to ALPERC.}
}
\references{
Arboretti, R., Ceccato, R., Pegoraro, L., Salmaso, L. () Active learning for noisy physical experiments with multiple responses
}
\author{
Luca Pegoraro
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{silhouette_plot}},\code{\link{ALPERC_2D_viz}}, \code{\link{ALPERC_2D_viz_multiv}}, \code{\link{ALPERC_3D_viz}}, \code{\link{ALPERC_3D_viz_multiv}}
}
\examples{
##Import required data
data(D_cand)
data(sigma_cand)
data(S)

##run the function
ALPERC(n_add=5,
       D_cand=D_cand,
       sigma_cand=sigma_cand,
       S=S,
       strategy="exploration",
       varimp_distance=T,
       n_clust=NULL,
       n_boot=500,
       alpha_rank=0.1,
       seed_rank=2105,
       paral=FALSE)

##run the function with parallel computation
library(doSNOW)

##initialize cluster
cl <- makeCluster(4)
registerDoSNOW(cl)

ALPERC(n_add=5,
       D_cand=D_cand,
       sigma_cand=sigma_cand,
       S=S,
       strategy="exploration",
       varimp_distance=T,
       n_clust=NULL,
       n_boot=500,
       alpha_rank=0.1,
       seed_rank=2105,
       paral=TRUE)

stopCluster(cl)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
