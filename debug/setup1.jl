using PyCall
emcee = pyimport("emcee")
using ForwardDiff
using WeightedEnsemble
using TestLandscapes
using BasicMD
using CovarianceMatrices

V = x -> SymmetricDoubleWell(x[1]);
cfg = ForwardDiff.GradientConfig(V, zeros(Float64, 1));
∇V! = (gradV, X) -> ForwardDiff.gradient!(gradV, V, X, cfg);

fB = X -> Float64(X[1] > b);

sampler = EM(∇V!, β, Δt);

function restartA!(state::BasicMD.EMState)
    if (state.x[1] > b)
        @. state.x = x0
        ∇V!(state.∇V, x0)
    end
    state
end
constraints = Constraints(restartA!, trivial_constraint!, nΔt_recycle, nΔt_coarse);
mutation_opts = MDOptions(n_iters=nΔt_coarse, n_save_iters=nΔt_coarse);

mutation! = X -> sample_trajectory!(X, sampler, constraints, options=mutation_opts);

B0, bin_id, rebin! = setup_Voronoi_bins(x_voronoi);
E0 = Dirac_to_Ensemble(x0, n_particles)
rebin!(E0, B0, 0);

function sample_statistics(fB_trajectory, n_particles, ΔT_recycle; n_burn_in_vals=[0])

    n_burn = length(n_burn_in_vals)
    fB_mean = zeros(n_burn)
    fB_var = zeros(n_burn)
    IAT_sokal = zeros(n_burn)
    IAT_nw = zeros(n_burn)

    for (j, n_burn_in) in enumerate(n_burn_in_vals)
        fB_mean[j] = mean(fB_trajectory[n_burn_in+1:end])
        fB_var[j] = var(fB_trajectory[n_burn_in+1:end])
        IAT_sokal[j] = emcee.autocorr.integrated_time(fB_trajectory[n_burn_in+1:end])[1]
        IAT_nw[j] = lrvar(ParzenKernel{NeweyWest}(), [fB_trajectory[n_burn_in+1:end];;], prewhite=true)[1, 1] / fB_var[j]
    end

    return (fB_mean, fB_var, IAT_sokal, IAT_nw)
end