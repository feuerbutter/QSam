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
### Generation of a target sample
To generate a sample, one only needs to supply the appropriate parameters in the script, *main.m*. One needs to adjust the various properties of the reference sample, such as, *kappa_w*, the composition percentage of the Wishart sample, and to specify, *pop*, the posterior for the target sample. The target sample points are then collected in the variable *prob_points_accepted* and the acceptance rate is given in thevariable in *ar*. 

One is advised to first begin with a relatively small *N_total*, e.g. 1e4 or 1e5, and play with the reference sample parameters to achieve a considerably efficient acceptance rate before generating a large target sample.

The sample produced is probabilities in the probability simplex, and it could be transformed into density matrices in the state space by
```
rho = prob2Rho(prob,POM)
```
### Verification of the produced sample
To verify that the target sample produced indeed follows the target distrbution, one could run the script, *main_verification.m*. It is only possible for low dimension where the size of the sample could be aproximated accurately with a practically sized uniform sample to calculate the posterior theoretically.
 
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
