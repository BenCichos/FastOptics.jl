struct Tilt{T} <: AbstractOpticalElement where {T <: AbstractFloat}
    alpha::T
    theta::T

    function Tilt(alpha::T, theta::T) where {T <: AbstractFloat}
        new{T}(alpha, theta)
    end
end

theta(t::Tilt) = t.theta
alpha(t::Tilt) = t.alpha
angle(t::Tilt) = (theta(t), alpha(t))

export Tilt, theta, alpha, angle