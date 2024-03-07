struct WEsampler{TM, TS, TR, TA} <: AbstractWEsampler
    mutation!::TM
    selection!::TS
    rebin!::TR
    analysis!::TA
end

struct DistributedWEsampler{TM, TS, TR, TA} <: AbstractWEsampler
    mutation::TM
    selection!::TS
    rebin!::TR
    analysis!::TA
end

"""
    WEsampler(mutation!, selection!, rebin!; analysis! = trivial_analysis!)

### Arguments
* `mutation!` - mutation step
* `selection!` - selection step
* `rebin!` - rebinning step
### Optional Arguments
* `analysis! = trivial_analysis!` - analysis step
"""
function WEsampler(mutation!::TM, selection!::TS, rebin!::TR; analysis!::TA = trivial_analysis!) where {TM<:Function, TS<:Function, TR<:Function, TA<:Function}
    return WEsampler(mutation!, selection!, rebin!, analysis!)
end

function DistributedWEsampler(mutation::TM, selection!::TS, rebin!::TR; analysis!::TA = trivial_analysis!) where {TM<:Function, TS<:Function, TR<:Function, TA<:Function}
    return DistributedWEsampler(mutation, selection!, rebin!, analysis!)
end
