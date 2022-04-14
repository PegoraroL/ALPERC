\name{ALPERC_2D_viz}
\alias{ALPERC_2D_viz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ALPERC_2D_viz
}
\description{
A function that produces a 2D visualization of points in D_add_j and LAMBDA (from an ALPERC element), together with the mean response (or uncertainty) as predicted by the model. It currently works only for hetGP and ranger models. For an effective visualization, the trained models (independent variables), D_add_j and LAMBDA should consider data rescaled to 0-1.
}
\usage{
ALPERC_2D_viz(
              D_add_j,
              LAMBDA,
              model,
              viz_factors,
              fixed_factors,
              fixed_factor_levels,
              what_plot
              )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D_add_j}{
A tibble, that is an output of ALPERC() function.
}
  \item{LAMBDA}{
A tibble, that is an output of ALPERC() function.
}
  \item{model}{
	A hetGP or ranger (trained via caret) model.
}
  \item{viz_factors}{
A character vector that includes the names of 2 factors that will be used in the plot.
}
  \item{fixed_factors}{
A character vector that includes the names of the factors that will be fixed. They are all the input variables except from those included in viz_factors.
}
  \item{fixed_factor_levels}{
A numeric vector that includes the chosen levels of the fixed_factors. The order must correspond to the order in fixed_factors.
}
  \item{what_plot}{
A character vector that is either "mean" or "uncertainty". In the first case, the mean response is plotted, in the latter case the uncertainty is plotted.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A 2D contour plot (blue-low, yellow-high) that shows how the mean response or uncertainty (from the model) varies with respect to viz_factors, when the fixed_factors are constrained. The symbols in the plot correspond to the configurations proposed in LAMBDA, and the different symbols indicate the different clusters. The color of the symbols corresponds to the position in the rank, from white (low position, small uncertainty) to blue (high position, high uncertainty). The points in D_add_j are indicated with a red number, that shows the number of replicates of that configuration (a configuration is considered as replicated if another point shares the same configuration of viz_factors, as fixed_factors are not considered).
}
\references{
Arboretti, R., Ceccato, R., Pegoraro, L., Salmaso, L. () Active learning for noisy physical experiments with multiple responses


Binois, M., Gramacy, R. B., & Ludkovski, M. (2018). Practical Heteroscedastic Gaussian Process Modeling for Large Simulation Experiments, Journal of Computational and Graphical Statistics, 27(4), 808-821, https://doi.org/10.1080/10618600.2018.1458625

Binois, M., & Gramacy, R. B. (2021). hetGP: Heteroskedastic Gaussian Process Modeling and Sequential Design in R. Journal of Statistical Software, 98(13), 1-44. https://doi.org/10.18637/jss.v098.i13

Wright MN, Ziegler A (2017). ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R. Journal of Statistical Software, 77(1), 1–17. doi: 10.18637/jss

Kuhn, M. (2008). Building Predictive Models in R Using the caret Package. Journal of Statistical Software, 28(5), 1 - 26. doi:http://dx.doi.org/10.18637/jss.v028.i05
}
\author{
Luca Pegoraro
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ALPERC}}, \code{\link{ALPERC_2D_viz_multiv}}, \code{\link{ALPERC_3D_viz}}, \code{\link{ALPERC_3D_viz_multiv}}
}
\examples{
##Import required data
data(D_cand)
data(sigma_cand)
data(S)

##import saved models
library(hetGP)

model<-readRDS(system.file("data", "Robot_arm.rds", package="ALPERC"))

##run the ALPERC function
ALPERC_results<-ALPERC(n_add=5,
                    D_cand=D_cand,
                    sigma_cand=sigma_cand,
                    S=S,
                    strategy="exploration",
                    varimp_distance=TRUE,
                    n_clust=NULL,
                    n_boot=500,
                    alpha_rank=0.1,
                    seed_rank=2105)

##obtain the 2D visualization
ALPERC_2D_viz(D_add_j=ALPERC_results$D_add_j,
              LAMBDA=ALPERC_results$LAMBDA,
              model=model,
              viz_factors=c("B", "C"),
              fixed_factors=c("A","D","E","F"),
              fixed_factor_levels=c(0.5,0.5,0.5,0.5),
              what_plot="mean")
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
