let 
    β = 10.0
    δ = 0.1
    fA(x) = Float64(-δ < x[1] < δ)
    
    V(x) = SymmetricDoubleWell(x[1])
    ∇V! = (gradV, X) -> ForwardDiff.gradient!(gradV, V, X)
    Δt = 1e-2  # time step
    sampler = MALA(V, ∇V!, β, Δt)

    nΔt = 10 # number of fine time steps per coarse WE step

    # define the mutation mapping
    opts = MDOptions(n_iters=nΔt)
    mutation! = x -> sample_trajectory!(x, sampler, options=opts)
    
    using WeightedEnsemble
    voronoi_pts = [[x_] for x_ in LinRange(-1.5, 1.5, 13)]
    B0, bin_id, rebin! = setup_Voronoi_bins(voronoi_pts)
    n_we_steps = 10^2
    n_particles = 10^3

    # the following is somwhat arbitrary, but it spreads particles through state space
    E0 = Ensemble([[x_] for x_ in LinRange(-1.0, 1.0, n_particles)])
    rebin!(E0, B0, 0)
    
    uni_selection! = (E, B, t) -> uniform_selection!(E, B, t)
    uni_we_sampler = WEsampler(mutation!, uni_selection!, rebin!)
    
    Random.seed!(1000)
    f_uni_trajectory = run_we_observables(E0, B0, uni_we_sampler, n_we_steps, (fA,))[:]
    p_est = mean(f_uni_trajectory[26:n_we_steps]);
    p_est ≈ 1.5255558001325998e-5
    
end