# common ensmble functions

function Base.length(E::TE) where {TE<:AbstractEnsembleType}
   return length(E.ξ)
end

function Base.isempty(E::TE) where {TE<:AbstractEnsembleType}
   return isempty(E.ξ)
end
