## Supplementary material for "Robustifying Bayesian nonparametric mixtures for count data" 

This repository contains 3 files:  

1) The current README file.
2) run.R:  The main R script implementing our model for the Okaloosa darter dataset (Dorazio et al., 2008)
3) RBNP4C_0.1.tar.gz: A binary R package file to be compiled locally before running run.R from R.
   
To compile RBNP4C_0.1.tar.gz simply write in any unix-like environment:
   
```
R CMD COMPILE RBNP4C_0.1.tar.gz 
```   

## References

Canale, A., Pruenster, I., (2016), Robustifying Bayesian nonparametric mixtures for count data, Biometrics, in press, doi:10.1111/biom.12538

Dorazio et al., (2008), "Modeling unobserved sources of heterogeneity in animal abundance using a Dirichlet	process prior", Biometrics 64, 635-644
