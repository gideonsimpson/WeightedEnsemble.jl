# Mutation

The mutation step is the user defined function corresponding to one step of a
(time-independent) Markov chain, ``X_{k+1}\sim K(X_k,\bullet)``.  This should be
coded as an in place transform, `mutation!`, such that one can execute
```
mutation!.(E.ξ);
```


## Molecular Dynamics Exmaple
The following example formats one time step of the Euler-Maruyama molecular
dynamics (MD) integrator in a way that it can be used with `WeightedEnsemble`.
This employs the [BasicMD](https://github.com/gideonsimpson/BasicMD.jl) package.  Here, we discretize the SDE,
```math
dX_t = -\nabla V(X_t)dt + \sqrt{2\beta^{-1}}dW_t,
```
with time step ``Δt``, and then take `n` steps of it to correspond to one step of the Markov chain:
```
sampler = EM(gradV!, β, Δt);
opts = MDOptions(n_iters = n);
mutation! = x -> sample_trajectory!(x, sampler, options = opts);
```
In a slight abuse of notation, the Markov chain, ``(X_k)`` that is used
within WE corresponds to the skeleton of the continuous in time process, ``(X_{k
n \Delta t })``.