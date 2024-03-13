# Resampling Algorithms

As it is neccessary to perform resampling, both for determining how many
offspring each bin should have, and, having determined that, resampling the
particles within each bin, we have provided several choices here.  These
methods, and their properties, can be found in Douc et al. (2005).  The provided
`selection!` functions default to using [`WeightedEnsemble.systematic`](@ref)
across bins and [`WeightedEnsemble.multinomial`](@ref) within bins.

All of these resampling methods take as their arguments, `n`, the number of
trials, and `ω`, a vector of probabilities.  What is returned is a vector,
`Nvals`, the same size as `ω`, containing the number allocated to each of the
states.  These will satisfy the properties that `sum(Nvals)=n`, and the mean
value of `Nvals[i]` (from repeated independent trials), is `n * ω[i]` for each `i`.

```@docs
WeightedEnsemble.residual
WeightedEnsemble.stratified
WeightedEnsemble.systematic
WeightedEnsemble.multinomial
```