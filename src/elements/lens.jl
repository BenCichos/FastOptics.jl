struct ThinLens{T} <: AbstractOpticalElement where {T <: AbstractFloat}
    f::T
    function ThinLens(f::T) where {T <: AbstractFloat}
        new{T}(f)
    end
end

focal_length(l::ThinLens) = l.f

export ThinLens, focal_length