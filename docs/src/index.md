# WeightedEnsemble.jl

```@contents
Pages = ["index.md"]
```


## Overview
Weighted Ensemble (WE) is a variance reduction strategy for computing rare
events.  One of the major advantages of WE over other variance reduction
strategies is that it makes rather minimal assumptions on the Markov underlying
process.  Indeed, there is no assumption of reversibility of the process, nor is
there any need to evaluate or otherwise estimate a likelihood.  It is sufficient
to be able to simulate the associated process.  A key application of this
approach is in the computation of mean first passage times (MFPTs) via the Hill
relation.

This implementation of WE makes use of a fixed number, ``N``, particles, ``\xi^{(i)}``, each carrying weight ``\omega{(i)}``.  At algorithmic time ``t=0``, we denote
```math
\mu_0^{(N)}(dx) = \sum_{i=1}^N \omega_0^{(i)}\delta_{\xi_0^{(i)}}(dx)
```
The weights are assumed to be positive and sum to one.  The particle ensemble is
then evolved in two steps

### Selection
During the selection step, particles are resampled from the empirical
distribution in such a way so as to maintain unbiasedness _and_ reduce variance
with respect to some quantity of interest.  The ensemble after selection is
represented by
```math
\hat{\mu}_0^{(N)}(dx) = \sum_{i=1}^N \hat{\omega}_0^{(i)}\delta_{\hat{\xi}_0^{(i)}}(dx)
```

### Mutation
During the mutation step, the particles evolve freely over one algorithmic time unit according to a user defined Markov kernel, $K$, with ``\xi_{1}^{(i)}\sim K(\hat{\xi}_0^{(0)}, dx)``, and ``\omega_{1}^{(i)}= \hat{\omega}_0^{(0)}``.  We then have the updated ensemble,
```math
\mu_1^{(N)}(dx) = \sum_{i=1}^N \omega_1^{(i)}\delta_{\xi_1^{(i)}}(dx)
```


### Unbiasedness 


## Caveats

## Collaborators
* D. Aristoff
* J. Copperman
* L. F. Doherty
* F. G. Jones
* R. J. Webber
* D. M. Zuckerman

## Acknowledgements
This work was supported in part by the US National Science Foundation Grants DMS-1818716 and DMS-2111278.

## References
 1.  Aristoff, D., Copperman, J., Simpson, G., Webber, R. J. & Zuckerman, D. M. Weighted ensemble: Recent mathematical developments. J. Chem. Phys. 158, 014108 (2023).
 2.  Russo, J. D. et al. WESTPA 2.0: High-Performance Upgrades for Weighted Ensemble Simulations and Analysis of Longer-Timescale Applications. J. Chem. Theory Comput. 18, 638–649 (2022).
 3.  Aristoff, D. An ergodic theorem for the weighted ensemble method. arXiv:1906.00856 (2021).
 4.  Webber, R. J., Aristoff, D., & Simpson, G. A splitting method to reduce MCMC variance. arXiv:2011.13899 (2020)
 5.  Aristoff, D. & Zuckerman, D. M. Optimizing Weighted Ensemble Sampling of Steady States. Multiscale Model. Simul. 18, 646–673 (2020).
 6.  Huber, G. A. & Kim, S. Weighted-ensemble Brownian dynamics simulations for protein association reactions. Biophysical Journal 70, 97–110 (1996).

 