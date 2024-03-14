# Data Structures for Weighted Ensemble

```@contents
Pages = ["struct1.md"]
```

## Ensembles
```@docs
Ensemble
```

The ensemble data structure holds the information about ``\{(\omega_t^{(i)}, \xi_t^{(i)}\}``, along with the post selection step ``\{(\hat{\omega}_t^{(i)}, \hat{\xi}_t^{(i)}\}``.  Fo convenience, it also carries information about the bin the particle currently resides in.  The structures `d` and `d̂` are optional for carrying any auxiliary information.

Several convenience constructors have been included.  
```@docs
Ensemble(X::TP, n_particles::TI) where {TP, TI<:Integer}
Ensemble(X::Vector{TP}) where {TP}
Ensemble(X::Vector{TP}, ω::Vector{TF}) where {TP, TF<:AbstractFloat}
```


As an example, if our problem is posed in ``\mathbb{R}``, an initial ensemble could be constructed:
```
x0 = [-1.0];
n_particles = 100;
E0 = Ensemble(x0, n_particles);
```
This initializes the bin id field, `b`, to zeros.  It will be
neccessary to either manually assign this field, or to call the user defined
`rebin!` function:
```
rebin!(E0, B0, 0);
```

## Bins
```@docs
Bins
```
In addition to defining the `Bins` data strucutre, it will also be essential to provide:
* `bin_id(x)` - This maps points in the state space to an integer indexing fo the bins.
* `rebin!(E, B, t)` - This determines the bin of each particle, records this in the ensemble `E.b` field, sums the weights of all particles within the bins, and updates this in the bins `B.ν` field.  This may be time dependent.

Several tools have been provided to simplify the construction of the bins.  

```@docs
WeightedEnsemble.Points_to_Bins
```
This is useful if the bins can be sensibily identified with a set of
points in the state space, as in the case of a Voronoi tesselation.  Indeed, as
Voronoi tesselations are often used, we have
```@docs
setup_Voronoi_bins
```
As an example, suppose our Voronoi bins are points in ``\mathbb{R}^2``, then the following code will produce all the needed functions:
```
voronoi_centers = [[-1.0, -2.0], [2.0, 1.0], [1.0, 0.]];
B0, bin_id, rebin! = setup_Voronoi_bins(voronoi_centers);
```

## Samplers
```@docs
WEsampler
```

