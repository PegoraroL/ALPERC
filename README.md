# ALPERC
This package provides the code implementing ALPERC: Active Learning for Physical Experiments based on nonparametric Ranking and Clustering. ALPERC can be used to drive batch sequential data acquisition when more than 2 responses are investigated in the same experiment and the objective of the study is to build models that maximize the predictive accuracy with respect to all the responses. The package also includes some functions for a 2D and 3D interactive visualization of the response surface (mean and uncertainty) developed by the models, and shows the experimental configurations selected by ALPERC.
For further information and examples, check package documentation and vignettes.

## Installation

devtools::install_github("PegoraroL/ALPERC")

## References

Arboretti, R., Ceccato, R., Pegoraro, L., Salmaso, L. (2022), Active learning for noisy physical experiments with more than two responses, Chemometrics and Intelligent Laboratory Systems, https://doi.org/10.1016/j.chemolab.2022.104595
