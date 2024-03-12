# Allocation


```@contents
Pages = ["allocation1.md"]
```

The allocation methods determine:
* How many particles we should attempt to have within each bin.  
* How many offspring each particle should have. 
In the methods included here,
one first allocates the number of particles to each bin, and then uses
multinomial resampling, within each bin, to determine how many offspring each
particle should have.  

## Bin Allocation Methods

```@docs
    minimal_bin_allocation!
```

```@docs
    uniform_bin_allocation!
```

```@docs
    optimal_bin_allocation!
```
This function is defined in the spirit of the optimal allocation analysis of Aristoff & Zuckerman (2020).

```@docs
    targeted_bin_allocation!
```

## Particle Allocation Methods
```@docs
    within_bin_allocation!
```

## Utility Functions
```@docs
    trivial_allocation!
```