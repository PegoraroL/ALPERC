# ALPERC
This package provides the code implementing ALPERC: Active Learning for Physical Experiments based on nonparametric Ranking and Clustering. ALPERC can be used to drive batch sequential data acquisition when more than 2 responses are investigated in the same experiment and the objective of the study is to build models that maximize the predictive accuracy with respect to all the responses. The package also includes some functions for a 2D and 3D interactive visualization of the response surface (mean and uncertainty) resulting from the models, and shows the experimental configurations selected by ALPERC.
For further information and examples, check package documentation and vignettes.

## Installation
Install from GitHub:
``` r
devtools::install_github("PegoraroL/ALPERC")
```

ALPERC: vignette
================


## **ALPERC()**

### **Description**

ALPERC is an Active Learning (AL) algorithm for noisy physical
experiments with multiple responses based on nonparametric ranking and
clustering.

For further information on the algorithm, please refer to the article
“Arboretti, R., Ceccato, R., Pegoraro, L., Salmaso, L. (2022), Active
learning for noisy physical experiments with more than two responses,
Chemometrics and Intelligent Laboratory Systems,
<https://doi.org/10.1016/j.chemolab.2022.104595>” and on the Help page
of the R package.

The functions have been implemented in the following environment, using
the following package versions:

``` r
#load the package
library(ALPERC)
```

    ## [1] "R version 4.1.1 (2021-08-10)"

    ##  ALPERC   hetGP forcats  tibble   tidyr   dplyr ggplot2 
    ## "0.1.0" "1.1.4" "0.5.1" "3.1.4" "1.1.4" "1.0.7" "3.3.5"

### **Usage**

ALPERC(n_add, D_cand, sigma_cand, S=NULL, strategy=“exploration”,
varimp_distance, n_clust=NULL, n_boot=1000, force_n\_boot=TRUE,
alpha_rank=0.1, seed_rank=2105, paral=FALSE)

### **Arguments**

| Argument        | Description                                                                                                                                                                                                                                                                                                                                                                                                                     |
|-----------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| n_add           | The size of the batch to add at the AL iteration. An integer.                                                                                                                                                                                                                                                                                                                                                                   |
| D_cand          | A data frame including the candidate experimental configurations for AL. The first column is a factor column that contains the identifier (ID) of the candidates, the following columns include the factor levels. The column headers (except from the first) are the names of the factors.                                                                                                                                     |
| sigma_cand      | A data frame including the predictive uncertainty of the candidate points in D_cand. The first column is a factor column that contains the identifier (ID) of the candidates, the following columns include the uncertainty of predictions of the candidate configurations for each response (>2). The column headers (except from the first) are the names of the responses.                                                   |
| S               | A data frame including the variable importance of each factor. The first column is a character column that contains the names of the factors (equal to the column headers in D_cand), the following columns include the variable importance of the factors corresponding to each response (>2). The column headers (except from the first) are the names of the responses. If variable importance is not available, set S=NULL. |
| strategy        | It is either “exploration” (only unique candidates) or “exploitation” (can include replicates).                                                                                                                                                                                                                                                                                                                                 |
| varimp_distance | It is either TRUE or FALSE. In the first case, the assignment of each candidate configuration to a cluster is adjusted by the importance of each dimension. If varimp_distance=FALSE, S is not used.                                                                                                                                                                                                                            |
| n_clust         | The number of clusters, it can either be NULL or an integer. If n_clust=NULL, the best number of clusters is chosen with the Silhouette method, otherwise it is set equal to the number indicated.                                                                                                                                                                                                                              |
| n_boot          | The number of bootstrap iterations for the permutations used for the nonparametric ranking procedure.                                                                                                                                                                                                                                                                                                                           |
| force_n\_boot   | Either TRUE or FALSE. If TRUE, it forces n_boot random permutations for computation of the nonparametric ranking. If FALSE, the function automatically evaluates whether n_boot is larger than the cardinality. If force_n\_boot=FALSE and n_boot \> cardinality, the function generates the exact permuted samples, speeding up the computation. Recommended only for 5 or more responses.                                     |
| alpha_rank      | The significance level for the computation of the nonparametric ranking procedure.                                                                                                                                                                                                                                                                                                                                              |
| seed_rank       | A seed that is used to initialize the permutations.                                                                                                                                                                                                                                                                                                                                                                             |
| paral           | Either TRUE or FALSE. It indicates whether parallel computation should be employed (default is FALSE). If paral=TRUE, a cluster must be initialized (e.g. via “doSNOW” library)                                                                                                                                                                                                                                                 |

### **Outputs**

The function provides a list with 8 elements:

| Element              | Description                                                                                                                                                       |
|----------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| $strategy            | Is the strategy chosen for ALPERC implementation.                                                                                                                 |
| $seed_rank           | Is the seed chosen to initialize the nonparametric ranking procedure.                                                                                             |
| $n_permut_comparisons | Is the number of permutations for the ranking procedure. If force_n\_boot=FALSE, it can be different from n_boot.                                                 |
| $varimp_distance     | Indicates whether the variable importance is used or not.                                                                                                         |
| $clust_choice        | A tibble that includes the number of clusters and corresponding Silhouette index.                                                                                 |
| $nclust_criterion    | Indicates the strategy used for the selection of the number of clusters.                                                                                          |
| $best_nclust         | The number of clusters considered.                                                                                                                                |
| $LAMBDA              | A tibble that includes the ID of each candidate, the factor levels combination, the position in the rank and the cluster to which each configuration is assigned. |
| $D_add_j             | A tibble that includes the n_add configurations from D_add that should be added according to ALPERC.                                                              |

### **Example**

Here we show an example of the application of ALPERC function. First, we
load the input data frames:

``` r
D_cand
```

    ## # A tibble: 50 x 7
    ##    rowid     A     B     C     D     E     F
    ##    <fct> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ##  1 1       0.6   0.2   0     0.2   0.2   0  
    ##  2 2       0.2   0.2   0.6   1     0.6   0.4
    ##  3 3       1     0.4   0.8   1     0.4   0.6
    ##  4 4       0     0.6   1     0.6   0.6   0.2
    ##  5 5       0.4   0.2   0     0.4   0.8   0.2
    ##  6 6       0.4   0     0.4   0.8   0     0.6
    ##  7 7       0     0.8   0     0     0.4   0.8
    ##  8 8       0     0.2   0.4   0.6   0.4   0  
    ##  9 9       0.4   0.6   0.2   0.8   1     0  
    ## 10 10      0     0.8   0.2   1     0.8   0.8
    ## # ... with 40 more rows

``` r
sigma_cand
```

    ## # A tibble: 50 x 8
    ##    rowid X6.d.Borehole X6.d.OTL.circuit X6.d.Piston X6.d.Piston.Mod X6.d.Robot.arm X6.d.Rosenbrock.Function
    ##    <fct>         <dbl>            <dbl>       <dbl>           <dbl>          <dbl>                    <dbl>
    ##  1 1            0.0245           0.0327      0.0837          0.0836         0.148                    0.0607
    ##  2 2            0.0185           0.0440      0.0718          0.0589         0.0834                   0.0989
    ##  3 3            0.0582           0.0177      0.0549          0.0563         0.0869                   0.122 
    ##  4 4            0.0175           0.0989      0.0676          0.0626         0.100                    0.132 
    ##  5 5            0.0218           0.0342      0.0696          0.0679         0.110                    0.0639
    ##  6 6            0.0353           0.0411      0.0928          0.0715         0.0894                   0.0595
    ##  7 7            0.0345           0.0968      0.0486          0.0627         0.137                    0.0832
    ##  8 8            0.0132           0.0642      0.0900          0.0717         0.0877                   0.0661
    ##  9 9            0.0268           0.0753      0.0459          0.0498         0.0926                   0.145 
    ## 10 10           0.0289           0.0423      0.0378          0.0398         0.107                    0.123 
    ## # ... with 40 more rows, and 1 more variable: X6.d.Wing.weight <dbl>

``` r
S
```

    ## # A tibble: 6 x 8
    ##   X     X6.d.Borehole X6.d.OTL.circuit X6.d.Piston X6.d.Piston.Mod X6.d.Robot.arm X6.d.Rosenbrock.Function
    ##   <chr>         <dbl>            <dbl>       <dbl>           <dbl>          <dbl>                    <dbl>
    ## 1 A           0.585            0.529       0.0239          0.0656          0.0110                  0.175  
    ## 2 B           0.00204          0.360       0.606           0.483           0.431                   0.191  
    ## 3 C           0.00187          0.0851      0.358           0.366           0.428                   0.198  
    ## 4 D           0.0409           0.0247      0.0202          0.162           0.149                   0.205  
    ## 5 E           0.0304           0.00387     0.00248         0.00128         0.231                   0.239  
    ## 6 F           0.478            0.00178     0.00507         0.00300         0.176                   0.00984
    ## # ... with 1 more variable: X6.d.Wing.weight <dbl>

Then, we use the function:

``` r
##Here we run the function with parallel computation. This step can be avoided by setting the argument paral=FALSE.
library(doSNOW)
```

``` r
##initialize cluster
cl <- makeCluster(4)
registerDoSNOW(cl)
                   
results<-ALPERC(n_add=5,
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
```

    ## [1] "Start NPC ranking"
    ## [1] "2022-07-04 10:44:10 CEST"
    ##|==================================================================================================================| 100%
    ## Time difference of 32.26095 secs
    ## [1] "End NPC ranking"

``` r
stopCluster(cl)
```

Finally, we print the results:

``` r
results
```

    ## $strategy
    ## [1] "exploration"
    ## 
    ## $seed_rank
    ## [1] 2105
    ## 
    ## $n_permut_comparisons
    ## [1] 500
    ## 
    ## $varimp_distance
    ## [1] TRUE
    ## 
    ## $clust_choice
    ## # A tibble: 25 x 2
    ##    clusters     y
    ##    <fct>    <dbl>
    ##  1 1        0    
    ##  2 2        0.149
    ##  3 3        0.118
    ##  4 4        0.129
    ##  5 5        0.143
    ##  6 6        0.149
    ##  7 7        0.158
    ##  8 8        0.161
    ##  9 9        0.152
    ## 10 10       0.157
    ## # ... with 15 more rows
    ## 
    ## $nclust_criterion
    ## [1] "silhouette"
    ## 
    ## $best_nclust
    ## [1] 18
    ## 
    ## $LAMBDA
    ## # A tibble: 50 x 9
    ##    rowid     A     B     C     D     E     F  rank cluster
    ##    <chr> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <int>   <int>
    ##  1 26      0.8   0.8   1     0     1     0.2    50      10
    ##  2 25      0     0.2   1     0.2   0.8   0.6    49       9
    ##  3 37      0.6   0.4   1     0     0.6   0.6    48      15
    ##  4 38      0.6   0     0     0.2   0.4   1      47      16
    ##  5 7       0     0.8   0     0     0.4   0.8    46       1
    ##  6 1       0.6   0.2   0     0.2   0.2   0      30       1
    ##  7 4       0     0.6   1     0.6   0.6   0.2    30       2
    ##  8 9       0.4   0.6   0.2   0.8   1     0      30       1
    ##  9 11      0.8   0     1     0.2   0.2   0.6    30       4
    ## 10 15      1     0.4   1     0.8   1     0.6    30       1
    ## # ... with 40 more rows
    ## 
    ## $D_add_j
    ## # A tibble: 5 x 7
    ##   rowid     A     B     C     D     E     F
    ##   <fct> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1 26      0.8   0.8     1   0     1     0.2
    ## 2 25      0     0.2     1   0.2   0.8   0.6
    ## 3 37      0.6   0.4     1   0     0.6   0.6
    ## 4 38      0.6   0       0   0.2   0.4   1  
    ## 5 7       0     0.8     0   0     0.4   0.8

## **silhouette_plot()**

A function that gets as input an “ALPERC” element (list) and provides a
Silhouette plot, indicating the best number of clusters.

### **Example**

``` r
silhouette_plot(results) 
```

![<https://imgur.com/Edn066h>](https://imgur.com/Edn066h.png)

## **ALPERC_2D_viz()**

A function that produces a 2D visualization of D_add_j and LAMBDA,
together with the mean response (or uncertainty) as predicted by the
model. It currently works only for hetGP and ranger models.

### **Usage**

ALPERC_2D_viz( D_add_j, LAMBDA, model_fnct, viz_factors, fixed_factors,
fixed_factor_levels, what_plot )

### **Arguments**

| Argument            | Description                                                                                                                                                |
|---------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------|
| D_add_j             | A tibble, that is the output of ALPERC() function.                                                                                                         |
| LAMBDA              | A tibble, that is the output of ALPERC() function.                                                                                                         |
| model_fnct          | A function that returns a list of two elements, where the first element indicates the predictive mean and the second the predictive uncertainty.           |
| viz_factors         | A character vector that includes the names of 2 factors that will be used in the plot.                                                                     |
| fixed_factors       | A character vector that includes the names of the factors that will be fixed. They are all the input variables except from those included in viz_factors.  |
| fixed_factor_levels | A numeric vector that includes the chosen levels of the fixed_factors. The order must correspond to the order in fixed_factors                             |
| what_plot           | A character vector that is either “mean” or “uncertainty”. In the first case, the mean response is plotted, in the latter case the uncertainty is plotted. |

### **Outputs**

A 2D contour plot (blue-low, yellow-high) that shows how the mean
response or uncertainty (from the model) varies with respect to
viz_factors, when the fixed_factors are constrained. The symbols in the
plot correspond to the configurations proposed in LAMBDA, and the
different symbols indicate the different clusters. The color of the
symbols corresponds to the position in the rank, from white (low
position, small uncertainty) to blue (high position, high uncertainty).
The points in D_add_j are indicated with a red number, that shows the
number of replicates of that configuration (a configuration is
considered as replicated if another point shares the same configuration
of viz_factors, as fixed_factors are not considered).

### **Example**

``` r
#here load the models with readRDS(). We load 7 models
```

``` r
#we consider hetGP model
model_fnct<-function(X){
    pred_m<-predict(x = X, object = model5)$mean
    pred_unc<-(sqrt(predict(x = X, object = model5)$nugs +
                        predict(x = X, object = model5)$sd2))
    return(list(pred_m,pred_unc))
}

plot2D_mean<-ALPERC_2D_viz(D_add_j=results$D_add_j,
                      LAMBDA=results$LAMBDA, 
                      model_fnct=model_fnct, 
                      viz_factors=c("B", "C"),
                      fixed_factors=c("A","D","E","F"), 
                      fixed_factor_levels=c(0.5,0.5,0.5,0.5),
                      what_plot="mean")
```
![https://imgur.com/fgSi82v](https://imgur.com/fgSi82v.png)

``` r
plot2D_uncert<-ALPERC_2D_viz(D_add_j=results$D_add_j,
                      LAMBDA=results$LAMBDA, 
                      model_fnct=model_fnct, 
                      viz_factors=c("B", "C"),
                      fixed_factors=c("A","D","E","F"), 
                      fixed_factor_levels=c(0.5,0.5,0.5,0.5),
                      what_plot="uncertainty")
```
![https://imgur.com/q7UKJ9I](https://imgur.com/q7UKJ9I.png)


## **ALPERC_2D_viz_multiv()**

A function equivalent to ALPERC_2D_viz(), but provides a visualization
of multiple responses in one plot. The contour plot is obtained as the
average of the mean responses (or uncertainty) from the models. For this
reason, raw data (responses) should be transformed to 0-1 before fitting
the models. It currently works only for hetGP and ranger models.

### **Usage**

ALPERC_2D_viz_multiv( D_add_j, LAMBDA, model_fnct, viz_factors,
fixed_factors, fixed_factor_levels, what_plot )

### **Arguments**

| Argument            | Description                                                                                                                                                                |
|---------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| D_add_j             | A tibble, that is the output of ALPERC() function.                                                                                                                         |
| LAMBDA              | A tibble, that is the output of ALPERC() function.                                                                                                                         |
| model_fnct          | A list of functions each returning a list of two elements, where the first element indicates the predictive mean and the second the predictive uncertainty.                |
| viz_factors         | A character vector that includes the names of 2 factors that will be used in the plot.                                                                                     |
| fixed_factors       | A character vector that includes the names of the factors that will be used fixed. They are all the input variables except from those included in viz_factors.             |
| fixed_factor_levels | A numeric vector that includes the chosen levels of the fixed_factors. The order must correspond to the order in fixed_factors                                             |
| what_plot           | A character vector that is either “mean” or “uncertainty”. In the first case, the average mean response is plotted, in the latter case the average uncertainty is plotted. |

### **Outputs**

A 2D contour plot (blue-low, yellow-high) that shows how the average
response or average uncertainty (over the responses as predicted by the
models) varies with respect to viz_factors, when the fixed_factors are
constrained. The symbols in the plot correspond to the configurations
proposed in LAMBDA, and the different symbols indicate the different
clusters. The color of the symbols corresponds to the position in the
rank, from white (low position, small uncertainty) to blue (high
position, high uncertainty).

### **Example**

``` r
#we consider hetGP models
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
```

``` r
plot2D_multiv_m<-ALPERC_2D_viz_multiv(D_add_j=results$D_add_j,
                      LAMBDA=results$LAMBDA, 
                      model_fnct=model_fnct_list, 
                      viz_factors=c("B", "C"),
                      fixed_factors=c("A","D","E","F"), 
                      fixed_factor_levels=c(0.5,0.5,0.5,0.5),
                      what_plot="mean")
```

![https://imgur.com/NsZr8dD](https://imgur.com/NsZr8dD.png)

``` r
plot2D_multiv_u<-ALPERC_2D_viz_multiv(D_add_j=results$D_add_j,
                      LAMBDA=results$LAMBDA, 
                      model_fnct=model_fnct_list, 
                      viz_factors=c("B", "C"),
                      fixed_factors=c("A","D","E","F"), 
                      fixed_factor_levels=c(0.5,0.5,0.5,0.5),
                      what_plot="uncertainty")
```
![https://imgur.com/Hfqm8XG](https://imgur.com/Hfqm8XG.png)

## **ALPERC_3D_viz()**

A function that provides a 3D visualization of D_add_j and LAMBDA (from
an ALPERC element), together with the mean response or uncertainty as
predicted by the model. It is equivalent to ALPERC_2D_viz(), but
produces an interactive surface plot in 3 dimensions. It currently works
only for hetGP and ranger models.

### **Usage**

ALPERC_3D_viz( D_add_j, LAMBDA, model_fnct, viz_factors, fixed_factors,
fixed_factor_levels, what_plot, n_pred_points )

### **Arguments**

| Argument            | Description                                                                                                                                                |
|---------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------|
| D_add_j             | A tibble, that is the output of ALPERC() function.                                                                                                         |
| LAMBDA              | A tibble, that is the output of ALPERC() function.                                                                                                         |
| model_fnct          | A function that returns a list of two elements, where the first element indicates the predictive mean and the second the predictive uncertainty.           |
| viz_factors         | A character vector that includes the names of 2 factors that will be used in the plot.                                                                     |
| fixed_factors       | A character vector that includes the names of the factors that will be fixed. They are all the input variables except from those included in viz_factors.  |
| fixed_factor_levels | A numeric vector that includes the chosen levels of the fixed_factors. The order must correspond to the order in fixed_factors                             |
| what_plot           | A character vector that is either “mean” or “uncertainty”. In the first case, the mean response is plotted, in the latter case the uncertainty is plotted. |
| n_pred_points       | An integer, that indicates the number of points (randomly sampled) from the surface, for which uncertainty of prediction is reported. It defaults to 0.    |

### **Outputs**

An interactive 3D surface plot that shows how the response (from the
model) or uncertainty varies with respect to viz_factors, when the
fixed_factors are constrained. In the bottom plane the candidate
configurations are shown with different colors, corresponding to the
different clusters. A number indicates the number of replicates (only
viz_factors are considered). Red points indicate the mean uncertainty of
predictions in some locations.

### **Example**

``` r
#we consider hetGP models
model_fnct<-function(X){
    pred_m<-predict(x = X, object = model5)$mean
    pred_unc<-(sqrt(predict(x = X, object = model5)$nugs +
                        predict(x = X, object = model5)$sd2))
    return(list(pred_m,pred_unc))
}



plot3D_m<-ALPERC_3D_viz(D_add_j=results$D_add_j,
                      LAMBDA=results$LAMBDA, 
                      model_fnct=model_fnct, 
                      viz_factors=c("B", "C"),
                      fixed_factors=c("A","D","E","F"), 
                      fixed_factor_levels=c(0.5,0.5,0.5,0.5),
                      what_plot="mean",
                      n_pred_points=100)

plot3D_m
```
![https://imgur.com/KNvS6EW](https://imgur.com/KNvS6EW.gif)

``` r
plot3D_u<-ALPERC_3D_viz(D_add_j=results$D_add_j,
                      LAMBDA=results$LAMBDA, 
                      model_fnct=model_fnct, 
                      viz_factors=c("B", "C"),
                      fixed_factors=c("A","D","E","F"), 
                      fixed_factor_levels=c(0.5,0.5,0.5,0.5),
                      what_plot="uncertainty",
                      n_pred_points=0)

plot3D_u
```
![https://imgur.com/V3opsrb](https://imgur.com/V3opsrb.gif)


## **ALPERC_3D_viz_multiv()**

A function equivalent to ALPERC_3D_viz(), but provides a visualization
of multiple responses in one plot. The surface is obtained as the
average of the responses (or uncertainty) from the models. For this
reason, raw data (responses) should be transformed to 0-1 before fitting
the models. It is equivalent to ALPERC_2D_viz_multiv(), but produces an
interactive surface plot in 3 dimensions. It currently works only for
hetGP and ranger models.

### **Usage**

ALPERC_3D_viz_multiv( D_add_j, LAMBDA, model_fnct, viz_factors,
fixed_factors, fixed_factor_levels, what_plot, n_pred_points )

### **Arguments**

| Argument            | Description                                                                                                                                                                |
|---------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| D_add_j             | A tibble, that is the output of ALPERC() function.                                                                                                                         |
| LAMBDA              | A tibble, that is the output of ALPERC() function.                                                                                                                         |
| model_fnct          | A list of functions each returning a list of two elements, where the first element indicates the predictive mean and the second the predictive uncertainty. Maximum 10.    |
| viz_factors         | A character vector that includes the names of 2 factors that will be used in the plot.                                                                                     |
| fixed_factors       | A character vector that includes the names of the factors that will be used fixed. They are all the input variables except from those included in viz_factors.             |
| fixed_factor_levels | A numeric vector that includes the chosen levels of the fixed_factors. The order must correspond to the order in fixed_factors                                             |
| what_plot           | A character vector that is either “mean” or “uncertainty”. In the first case, the average mean response is plotted, in the latter case the average uncertainty is plotted. |
| n_pred_points       | An integer, that indicates the number of points (randomly sampled) from the surface, for which uncertainty of the prediction is reported. It defaults to 0.                |

### **Outputs**

An interactive 3D surface plot that shows how the mean response (average
from the models) or the mean uncertainty (average from the models)
varies with respect to viz_factors, when the fixed_factors are
constrained. In the bottom plane the candidate configurations are shown
with different colors, corresponding to the different clusters. A number
indicates the number of replicates (only viz_factors are considered).
Red points indicate the mean uncertainty of predictions (calculated as
the mean uncertainty from all models in the same location) in some
locations.

### **Example**

``` r
plot3D_multiv_m<-ALPERC_3D_viz_multiv(D_add_j=results$D_add_j,
                      LAMBDA=results$LAMBDA, 
                      model_fnct=model_fnct_list, 
                      viz_factors=c("B", "C"),
                      fixed_factors=c("A","E","D","F"), 
                      fixed_factor_levels=c(0.5,0.5,0.5,0.5),
                      what_plot="mean",
                      n_pred_points=100)

plot3D_multiv_m
```
![https://imgur.com/pEqWUax](https://imgur.com/pEqWUax.gif)

``` r
plot3D_multiv_u<-ALPERC_3D_viz_multiv(D_add_j=results$D_add_j,
                      LAMBDA=results$LAMBDA, 
                      model_fnct=model_fnct_list, 
                      viz_factors=c("B", "C"),
                      fixed_factors=c("A","E","D","F"), 
                      fixed_factor_levels=c(0.5,0.5,0.5,0.5),
                      what_plot="uncertainty",
                      n_pred_points=0)

plot3D_multiv_u
```
![https://imgur.com/rgByWKm](https://imgur.com/rgByWKm.gif)

## References

Arboretti, R., Ceccato, R., Pegoraro, L., Salmaso, L. (2022), Active learning for noisy physical experiments with more than two responses, Chemometrics and Intelligent Laboratory Systems, https://doi.org/10.1016/j.chemolab.2022.104595
