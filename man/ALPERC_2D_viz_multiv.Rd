\name{ALPERC_2D_viz_multiv}
\alias{ALPERC_2D_viz_multiv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ALPERC_2D_viz_multiv
}
\description{
A function equivalent to ALPERC_2D_viz(), but provides a visualization of multiple responses in one plot. The contour plot is obtained as the average of the mean responses (or uncertainty) from the models. For this reason, raw data (responses) should be transformed to 0-1 before fitting the models. For an effective visualization, the trained models (independent variables), D_add_j and LAMBDA should consider data rescaled to 0-1.
}
\usage{
ALPERC_2D_viz_multiv(
                    D_add_j,
                    LAMBDA,
                    model_fnct,
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
  \item{model_fnct}{
		A list of functions each returning a list of two elements, where the first element indicates the predictive mean and the second the predictive uncertainty.
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
A 2D contour plot (blue-low, yellow-high) that shows how the average response or average uncertainty (over the responses as predicted by the models) varies with respect to viz_factors, when the fixed_factors are constrained. The symbols in the plot correspond to the configurations proposed in LAMBDA, and the different symbols indicate the different clusters. The color of the symbols corresponds to the position in the rank, from white (low position, small uncertainty) to blue (high position, high uncertainty).
}
\references{
Arboretti, R., Ceccato, R., Pegoraro, L., Salmaso, L. (2022), Active learning for noisy physical experiments with more than two responses, Chemometrics and Intelligent Laboratory Systems, https://doi.org/10.1016/j.chemolab.2022.104595

Binois, M., Gramacy, R. B., & Ludkovski, M. (2018). Practical Heteroscedastic Gaussian Process Modeling for Large Simulation Experiments, Journal of Computational and Graphical Statistics, 27(4), 808-821, https://doi.org/10.1080/10618600.2018.1458625

Binois, M., & Gramacy, R. B. (2021). hetGP: Heteroskedastic Gaussian Process Modeling and Sequential Design in R. Journal of Statistical Software, 98(13), 1-44. https://doi.org/10.18637/jss.v098.i13

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
\code{\link{ALPERC}}, \code{\link{ALPERC_2D_viz}}, \code{\link{ALPERC_3D_viz}}, \code{\link{ALPERC_3D_viz_multiv}}
}
\examples{
##Import required data
data("D_cand")
data("sigma_cand")
data("S")

##import saved models
library(hetGP)

model1<-readRDS(system.file("data", "Borehole.rds", package="ALPERC"))
model2<-readRDS(system.file("data", "OTL_circuit.rds", package="ALPERC"))
model3<-readRDS(system.file("data", "Piston.rds", package="ALPERC"))
model4<-readRDS(system.file("data", "Piston_mod.rds", package="ALPERC"))
model5<-readRDS(system.file("data", "Robot_arm.rds", package="ALPERC"))
model6<-readRDS(system.file("data", "Rosenbrock_Function.rds", package="ALPERC"))
model7<-readRDS(system.file("data", "Wing_weight.rds", package="ALPERC"))


##define function that returns predictive mean and uncertainty given the model
model_fnct1<-function(X){
  pred_m<-predict(x = X, object = model1)$mean
  pred_unc<-(sqrt(predict(x = X, object = model1)$nugs +
                       predict(x = X, object = model1)$sd2))
  return(list(pred_m,pred_unc))
}
model_fnct2<-function(X){
  pred_m<-predict(x = X, object = model2)$mean
  pred_unc<-(sqrt(predict(x = X, object = model2)$nugs +
                       predict(x = X, object = model3)$sd2))
  return(list(pred_m,pred_unc))
}
model_fnct3<-function(X){
  pred_m<-predict(x = X, object = model3)$mean
  pred_unc<-(sqrt(predict(x = X, object = model3)$nugs +
                       predict(x = X, object = model3)$sd2))
  return(list(pred_m,pred_unc))
}
model_fnct4<-function(X){
  pred_m<-predict(x = X, object = model4)$mean
  pred_unc<-(sqrt(predict(x = X, object = model4)$nugs +
                       predict(x = X, object = model4)$sd2))
  return(list(pred_m,pred_unc))
}
model_fnct5<-function(X){
  pred_m<-predict(x = X, object = model5)$mean
  pred_unc<-(sqrt(predict(x = X, object = model5)$nugs +
                       predict(x = X, object = model5)$sd2))
  return(list(pred_m,pred_unc))
}
model_fnct6<-function(X){
  pred_m<-predict(x = X, object = model6)$mean
  pred_unc<-(sqrt(predict(x = X, object = model6)$nugs +
                       predict(x = X, object = model6)$sd2))
  return(list(pred_m,pred_unc))
}
model_fnct7<-function(X){
  pred_m<-predict(x = X, object = model7)$mean
  pred_unc<-(sqrt(predict(x = X, object = model7)$nugs +
                       predict(x = X, object = model7)$sd2))
  return(list(pred_m,pred_unc))
}


model_fnct_list<-list(model_fnct1,
                      model_fnct2,
                      model_fnct3,
                      model_fnct4,
                      model_fnct5,
                      model_fnct6,
                      model_fnct7
                      )

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
ALPERC_2D_viz_multiv(D_add_j=ALPERC_results$D_add_j,
                    LAMBDA=ALPERC_results$LAMBDA,
                    model_fnct=model_fnct_list,
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
