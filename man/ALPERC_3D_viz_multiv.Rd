\name{ALPERC_3D_viz_multiv}
\alias{ALPERC_3D_viz_multiv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ALPERC_3D_viz_multiv
}
\description{
A function equivalent to ALPERC_3D_viz(), but provides a visualization of multiple responses in one plot. The surface is obtained as the average of the responses (or uncertainty) from the models. For this reason, raw data (responses) should be transformed to 0-1 before fitting the models.  It is equivalent to ALPERC_2D_viz_multiv(), but produces an interactive surface plot in 3 dimensions. It currently works only for hetGP and ranger models. For an effective visualization, the trained models (independent variables), D_add_j and LAMBDA should consider data rescaled to 0-1.
}
\usage{
ALPERC_3D_viz_multiv(
                    D_add_j,
                    LAMBDA,
                    model,
                    viz_factors,
                    fixed_factors,
                    fixed_factor_levels,
                    what_plot,
                    n_pred_points
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
	A list of hetGP or ranger (trained via caret) models. Maximum 10.
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
  \item{n_pred_points}{
An integer, that indicates the number of points (randomly sampled) from the surface, for which uncertainty of prediction is reported. Default is 0.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An interactive 3D surface plot that shows how the mean response (average from the models) or the mean uncertainty (average from the models) varies with respect to viz_factors, when the fixed_factors are constrained. In the bottom plane the candidate configurations are shown with different colors, corresponding to the different clusters. A number indicates the number of replicates (only viz_factors are considered). Red points indicate the mean uncertainty of predictions (calculated as the mean uncertainty from all models in the same location) in some locations.
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ALPERC}}, \code{\link{ALPERC_2D_viz}}, \code{\link{ALPERC_2D_viz_multiv}}, \code{\link{ALPERC_3D_viz}}
}
\examples{
##Import required data
data(D_cand)
data(sigma_cand)
data(S)

##import saved models
model1<-readRDS("model1.rds")
model2<-readRDS("model2.rds")
model3<-readRDS("model3.rds")
model4<-readRDS("model4.rds")

model_list<-list(model1,
                 model2,
                 model3,
                 model4)


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

##obtain the interactive 3D visualization
ALPERC_3D_viz_multiv(D_add_j=ALPERC_results$D_add_j,
                    LAMBDA=ALPERC_results$LAMBDA,
                    model=model_list,
                    viz_factors=c("B", "C"),
                    fixed_factors=c("A","D","E","F"),
                    fixed_factor_levels=c(0.5,0.5,0.5,0.5),
                    what_plot="uncertainty",
                    n_pred_points=0)
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
