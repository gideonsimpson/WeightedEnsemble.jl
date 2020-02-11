"""
`Residual`: perform residual sampling

### Arguments
`n` - number of trials
`ω` - probabilities
"""
function Residual(n,ω)

    m = length(ω);

    R = sum(floor.(Int, n * ω));

    if (R < n)
        ω̄ = (n * ω .- floor.(n * ω)) / (n - R);
        N̄vals = rand(Multinomial(n-R, ω̄));
        Nvals = floor.(Int, n * ω) .+ N̄vals;
    else
        Nvals = floor.(Int, n * ω);
    end

    return Nvals

end

"""
`Stratified`: perform stratified sampling

### Arguments
`n` - number of trials
`ω` - probabilities
"""
function Stratified(n,ω)
    U = range(0,stop=n-1)/n .+ rand(n)/n;

    Nvals = counts(quantile.(Categorical(ω), U), 1:length(ω));

    return Nvals
end

"""
`Systematic`: perform systematic sampling

### Arguments
`n` - number of trials
`ω` - probabilities
"""
function Systematic(n,ω)

    U = range(0,stop=n-1)/n .+ rand()/n;

    Nvals = counts(quantile.(Categorical(ω), U), 1:length(ω));

    return Nvals
end
