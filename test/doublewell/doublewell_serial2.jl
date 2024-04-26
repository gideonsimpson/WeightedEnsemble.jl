let
    #=
    WE estimation of the probability of for a diffusion with X(0) = a satisfying
    X(T) ∈ (b, ∞) for the double well potential V(x) = (x²-1)².

    Test with optimal allocation.
    =#

    # include("doublewell_setup.jl");

    V(x) = SymmetricDoubleWell(x[1])
    ∇V! = (gradV, X) -> ForwardDiff.gradient!(gradV, V, X)

    a = -1.0   # starting point
    x₀ = [a]
    b = 0.5    # target set, (b, ∞)
    Δt = 1e-4  # time step
    β = 5.0    # inverse temperature
    T = 1.0    # terminal time
    nΔt = ceil(Int, T / Δt)

    f = x -> Int(x[1] > b) # define observable

    # define sampler
    sampler = MALA(V, ∇V!, β, Δt)


    # number of coarse steps in WE
    n_we_steps = 10;
    # number of time steps during mutation step
    nΔt_coarse = ceil(Int, nΔt / n_we_steps);
    # number of samples in coarse matrix
    n_samples_per_bin = 10^2;
    # ensemble size
    n_particles = 10^2;

    # define bin structure
    voronoi_pts = [[x] for x in LinRange(a - 0.1, b + 0.1, 21)];
    B₀, bin_id, rebin! = setup_Voronoi_bins(voronoi_pts);

    # define the mutation mapping
    opts = MDOptions(n_iters=nΔt_coarse, n_save_iters=nΔt_coarse)
    mutation! = x -> sample_trajectory!(x, sampler, options=opts);

    # construct coarse model matrix
    Random.seed!(100);
    x0_vals = copy(voronoi_pts);
    n_bins = length(B₀);
    K̃ = WeightedEnsemble.build_coarse_transition_matrix(mutation!, bin_id, x0_vals, n_bins, n_samples_per_bin);

    # define coarse observable as a bin function
    f̃ = f.(voronoi_pts);
    _, v²_vectors = WeightedEnsemble.build_coarse_vectors(n_we_steps, K̃, float.(f̃));
    v² = (x, t) -> v²_vectors[t+1][bin_id(x)]
    # define selection function
    selection! = (E, B, t) -> optimal_selection!(E, B, v², t)

    # set up sampler - trivial analysis step
    we_sampler = WEsampler(mutation!, selection!, rebin!);
    # set up ensemble
    E₀ = Ensemble(x₀, n_particles);
    rebin!(E₀, B₀, 0);
    # run
    E = deepcopy(E₀);
    B = deepcopy(B₀);
    Random.seed!(200)
    run_we!(E, B, we_sampler, n_we_steps);
    p_est = f.(E.ξ) ⋅ E.ω;
    p_est ≈ 0.0012233542285717757
end