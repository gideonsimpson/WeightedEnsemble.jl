# Weighted Ensemble Methods

```@contents
Pages = ["we1.md"]
```

## Serial Methods
```@docs
run_we
run_we_observables
run_we!
```

## Multithreated Methods
These methods make use of multithreading in the mutation step.

```@docs
trun_we
trun_we_observables
trun_we!
```
## Distributed Methods
These methods make use of distributed computation in the mutation step via `pmap`.

```@docs
prun_we
prun_we_observables
prun_we!
```
