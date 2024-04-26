let

    β = 10.0
    x0 = -1.0
    a = -2
    b = 0.5

    fB = X -> Float64(X[1] > b) # define observable

    V(x) = SymmetricDoubleWell(x[1])
    ∇V! = (gradV, X) -> ForwardDiff.gradient!(gradV, V, X)
    Δt = 1e-2  # time step
    Δt_recycle = 1e-2
    nΔt_recycle = Int(Δt_recycle / Δt) # number of time steps before applying recycler
    nΔt_coarse = 1 * nΔt_recycle # number of time steps in a coarse step

    sampler = EM(∇V!, β, Δt)

    mutation_opts = MDOptions(n_iters=nΔt_coarse, n_save_iters=nΔt_coarse)

    function restartA!(state::BasicMD.EMState)
        if (state.x[1] > b)
            @. state.x = [x0]
            ∇V!(state.∇V, [x0])
        end
        state
    end

    constraints = Constraints(restartA!, trivial_constraint!, nΔt_recycle, nΔt_coarse)

    mutation! = x -> sample_trajectory!(x, sampler, constraints, options=mutation_opts)

    Δb = 0.2
    x_voronoi = [[x - 0.5 * Δb] for x in -1.5:Δb:0.5+Δb]
    n_bins = length(x_voronoi)
    B0, bin_id, rebin! = setup_Voronoi_bins(x_voronoi)

    n_particles = 10^3
    E0 = Ensemble([[x_] for x_ in LinRange(-1.5, 0.5, n_particles)])
    rebin!(E0, B0, 0)

    uni_we_sampler = WEsampler(mutation!, uniform_selection!, rebin!)

    n_we_steps = 10^3

    Random.seed!(500);
    f_uni_trajectory = run_we_observables(E0, B0, uni_we_sampler, n_we_steps, (fB,))[:];
    p_est = mean(f_uni_trajectory[n_we_steps÷2+1:end]);
    p_est ≈ 3.9583318379425827e-7
end