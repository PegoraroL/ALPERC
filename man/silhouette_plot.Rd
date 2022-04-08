\name{silhouette_plot}
\alias{silhouette_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{silhouette_plot
}
\description{A function that gets as input the results from "ALPERC()" (list) and provides a Silhouette plot, indicating the best number of clusters.
}
\usage{
silhouette_plot(ALPERC_results)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ALPERC_results}{The list obtained as an output from ALPERC() function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A Silhouette plot with a vertical dashed line that indicates the best number of clusters according to the Silhouette index.
}
\references{
Arboretti, R., Ceccato, R., Pegoraro, L., Salmaso, L. () Active learning for noisy physical experiments with multiple responses

Rousseeuw, P. J. (1987). Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. Journal of computational and applied mathematics, 20, 53-65.
}
\author{
Luca Pegoraro
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ALPERC}}
}
\examples{
##Import required data
data(D_cand)
data(sigma_cand)
data(S)

##run the ALPERC function
ALPERC_results<-ALPERC(n_add=5,
                    D_cand=D_cand,
                    sigma_cand=sigma_cand,
                   S=S,
                    strategy="exploration",
                    varimp_distance=T,
                    n_clust=NULL,
                    n_boot=500,
                    alpha_rank=0.1,
                    seed_rank=2105)

##obtain the silhouette plot
silhouette_plot(ALPERC_results)

## The function is currently defined as
function (ALPERC_results)
{
    ALPERC_results$clust_choice \%>\% ggplot(aes(x = clusters,
        y = y, group = 1)) + geom_point() + geom_line() + geom_vline(aes(xintercept = ALPERC_results$clust_choice \%>\%
        filter(y == max(y)) \%>\% pull(clusters)), linetype = "dashed") +
        ylab("Average silhouette width") + xlab("Number of clusters") +
        theme_bw()
  }
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
