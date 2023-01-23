# One-Dimensional Functions

function frect end
function fsinc end
function ftriangle end
function fgaussian end

frect(x::T; w::T=1.0) where {T <: AbstractFloat} = abs(x) < w/2 ? 1 : 0
fsinc(x::T; w::T=1.0) where {T <: AbstractFloat} = sinc(x/(2*w))
ftriangle(x::T; w::T=1.0) where {T <: AbstractFloat} = (abs(x/(2*w)) > 1) ? 0 : 1 - abs(x/(2*w))
fgaussian(x::T; w::T=1.0) where {T <: AbstractFloat} = exp(-pi*x^2/(2*w))

export frect, fsinc, ftriangle, fgaussian

# Two-Dimensional Functions

function frect2d end
function fsinc2d end
function fgaussian2d end
function fcirc end
function fellispe end

frect2d(x::T, y::T; w::T=1.0, h::T=1.0) where {T <: AbstractFloat} = frect(x; w=w) * frect(y; w=h)
fsinc2d(x::T, y::T, w::T=1.0, h::T=1.0) where {T <: AbstractFloat} = fsinc(x; w=w) * fsinc(y; w=h)
fcircle(x::T, y::T; radius::T=1.0) where {T <: AbstractFloat} = (x^2 + y^2) < radius^2
fellipse(x::T, y::T; size::NTuple{2, T}) where {T <: AbstractFloat} = (( ( x/size[1] )^2 + ( y/size[2] )^2 ) < 1) ? 1 : 0
fellipse(x::T, y::T, radius_x::T, radius_y::T) where {T <: AbstractFloat} = (( ( x/radius_x )^2 + ( y/radius_y )^2 ) < 1) ? 1 : 0
fgaussian2d(x::T, y::T; size::NTuple{2, T}) where {T <: AbstractFloat} = fgaussian(x, w=size[1]) * fgaussian(y, w=size[2])
fgaussian2d(x::T, y::T, w::T, h::T) where {T <: AbstractFloat} = fgaussian(x, w=w) * fgaussian(y, w=h)


export frect2d, fsinc2d, fcirc, fellipse, fgaussian, fgaussian2d