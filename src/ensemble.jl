# common ensmble functions

function Base.length(E::TE) where {TE<:AbstractEnsemble}
   return length(E.ξ)
end

function Base.isempty(E::TE) where {TE<:AbstractEnsemble}
   return isempty(E.ξ)
end
