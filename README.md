# QSam
Generation of samples of quantum states according to target posterior according to the method outlined in 
*Uncorrelated problem-specific samples of quantum states from zero-mean Wishart distributions* by *R. Han et al.* 
([arXiv:2106.08533](https://arxiv.org/abs/2106.08533)). 

**bold**
*italic*

To include a quote:
> This is a quote

To include codes:
```
some code
```

To create a link:
This is from [name of link](https://address.of.the.site)

To make a list:
- sth
- sth
or
1. sth
2. sth

## Running of Codes
To generate a sample, one only needs to supply the appropriate parameters in the script, *main.m*, such as the posterior.

## Utility Functions
- genRefSam : used to generate reference sample

## Conventions
- When a sample is described as uniform, we mean that it is uniform in the sense of Hilbert-Schmidt distance.

## change of notations
- [x] d -> m
- [x] etaw -> kappa_w
- [x] etas -> kappa_s
- [x] etab -> kappa
- [x] CVM -> Sigma
- [x] A -> Psi
- [x] Ncolw -> n_w
- [x] Ncols -> n_s
- [x] rhop -> rho_peak
- [x] kra -> A 
- [x] corp -> prob_points
- [x] minlogx -> min_log_lambda
- [x] maxlogx -> max_log_lambda
- [x] dellogx -> log_lambda_step

## qMLE
https://github.com/qMLE/qMLE
the superfast MLE for non-physical posterior is from the paper, 
J. Shang, Z. Zhang, and H. K. Ng, "Superfast maximum likelihood reconstruction for quantum tomography," Phys. Rev. A 95, 062336 (2017).
It contains 
- qmt.m
- qse_apg.m
- qse_cgls.m
- proj_spectrahedron.m
