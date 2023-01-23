struct RectangleAperture{T}  <: FunctionalAperture where {T <: AbstractFloat}
    fn::Function
    size::NTuple{2, T}    

    function RectangleAperture(size::NTuple{2, T}) where {T <: AbstractFloat}
        any(size .<= 0) && throw(ArgumentError("The size of the aperture has to be positive."))
        fn(field::LightField) = rectangle_aperture(field, size) 
        new{T}(fn, size)
    end
end

RectangleAperture(width::T, height::T) where {T <: AbstractFloat} = RectangleAperture( (width, height) )
SquareAperture(size::T) where {T <: AbstractFloat} = RectangleAperture((size, size))

width(r::RectangleAperture) = r.size[1]
height(r::RectangleAperture) = r.size[2]
size(r::RectangleAperture) = r.size

function rectangle_aperture(field::LightField, size::NTuple{2, T}) where {T <: AbstractFloat}
    x_range, y_range = range(field)
    frect2d.(x_range, y_range', w=size[1], h=size[2])
end

export RectangleAperture, SquareAperture, width, height, size
export rectangle_aperture