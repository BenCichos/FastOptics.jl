struct CircleAperture{T} <: FunctionalAperture where {T <: AbstractFloat}
    fn::Function
    radius::T
    
    function CircleAperture(radius::T) where {T <: AbstractFloat}
        radius > 0 || throw(ArgumentError("The radius of the aperture has to be larger than zero"))
        fn(field::LightField) = circle_aperture(field, radius)
        new{T}(fn, radius)
    end
end


function circle_aperture(field::LightField, radius::T) where {T <: AbstractFloat}
    x_range, y_range = range(field)
    fcircle.(x_range, y_range', radius=radius)
end

struct EllipseAperture{T} <: FunctionalAperture where {T <: AbstractFloat}
    fn::Function
    size::NTuple{2, T}
    
    function EllipseAperture(size::NTuple{2, T}) where {T <: AbstractFloat}
        fn(field::LightField) = ellipse_aperture(field, size)
        new{T}(fn, size)
    end
end

EllipseAperture(width::T, height::T) where {T <: AbstractFloat} = EllipseAperture((width, height))

function ellipse_aperture(field::LightField, size::NTuple{2, T}) where {T <: AbstractFloat}
    x_range, y_range = range(field)
    fellipse.(x_range, y_range', size=size)
end

radius(ca::CircleAperture) = ca.radius
size(ea::EllipseAperture) = ea.size
width(ea::EllipseAperture) = ea.size[1]
height(ea::EllipseAperture) = ea.size[2]

export CircleAperture, EllipseAperture
export radius, size, width, height