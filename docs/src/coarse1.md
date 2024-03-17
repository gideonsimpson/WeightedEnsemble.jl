# Coarse Models

```@contents
Pages = ["coarse1.md"]
Depth = 2
```



Several tools are included for building the variance and discrepancy functions,
discussed in [aristoff_analysis_2018, aristoff_optimizing_2020,
webber_splitting_2020, aristoff_weighted_2023](@cite).  These methods require,
first, the construction of a coarse scale transition matrix ``\tilde{k}``.  This
typically corresponds to an approximation of the Markov kernel, ``K``, on the
user defined bins,
```math
\tilde{K}_{pq} \approx \mathbb{P}(p\to q),
```
Having, computed this matrix, we next solve for the coarse scale discrepancy and one step
mutation variance functions:
```math
\begin{gather}
(I - \tilde{K})\tilde{h} = \tilde{f} - \tilde{\mu}(\tilde{f})\\
\tilde{v}^2(p) = \mathrm{Var}_{\tilde{K}_{p,\bullet}}(\tilde{h})
\end{gather}
```

## Serial Methods
```@docs
WeightedEnsemble.build_coarse_transition_matrix
WeightedEnsemble.build_coarse_poisson
WeightedEnsemble.build_coarse_vectors
```
While `build_coarse_poisson` is appropriate when using WE with steady state problems, `build_coarse_vectors` is what should be invoked for finite time horizon problems; see [aristoff_analysis_2018](@cite).


## Multithreaded Methods
TBW

## Distributed Parellel Methods
TBW