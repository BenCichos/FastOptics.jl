struct CoherentCirclePupil{T} <: AbstractCoherentPupil where {T <: AbstractFloat}
    fn::Function
    radius::T

    function CoherentCirclePupil(radius::T, distance::T) where {T <: AbstractFloat}
        fn(field::LightField) = circle_pupil(field, radius, distance)
        new{T}(fn, radius)
    end
end


function circle_pupil(field::LightField, radius::NTuple{2, T}, distance::T=distance) where {T <: AbstractFloat}
    lambda = wavelength(field)
    fc = radius / (lambda * distance)

    fx_range, fy_range = fourier_range(field)
    H = fcircle.(fx_range, fy_range', fc)

    field_shifted = fftshift(matrix(field))
    fft!(field_shifted)
    field_shifted .*= H
    ifft!(field_shifted)
    ifftshift!(matrix(field), field_shifted)
end

struct CoherentEllipsePupil{T} <: AbstractCoherentPupil where {T <: AbstractFloat}
    fn::Function
    size::NTuple{2, T}

    function CoherentEllipsePupil(size::NTuple{2, T}, distance::T) where {T <: AbstractFloat}
        fn(field::LightField) = ellipse_pupil(field, size, distance)
        new{T}(fn, size)
    end 
end

function ellipse_pupil(field::LightField, size::NTuple{2, T}, distance::T) where {T <: AbstractFloat}
    lambda = wavelength(field)
    fc_x = size[1] / (lambda * distance)
    fc_y = size[2] / (lambda * distance)

    fx_range, fy_range = fourier_range(field)
    H = fellipse.(fx_range, fy_range', fc_x, fc_y)

    field_shifted = fftshift(matrix(field))
    fft!(field_shifted)
    field_shifted .*= H
    ifft!(field_shifted)
    ifftshift!(matrix(field), field_shifted)
end

struct CoherentRectanglePupil{T} <: AbstractCoherentPupil where {T <: AbstractFloat}
    fn::Function
    size::NTuple{2,T}

    function CoherentRectanglePupil(size::NTuple{2, T}, distance::T)  where {T <: AbstractFloat}
        fn(field::LightField) = rectangle_pupil(field, size, distance)
        new{T}(fn, size)
    end
end

function rectangle_pupil(field::LightField, size::NTuple{2, T}, distance::T) where {T <: AbstractFloat}
    x_size, y_size = fieldsize(field)
    x_real, y_real = fieldsize(field)
    lambda = wavelength(field)
    fc_x = size[1] / (lambda * distance)
    fc_y = size[2] / (lambda * distance)

    fx = fftfreq(x_size, x_size / x_real)
    fy = fftfreq(y_size, y_size / y_real)
    H  = frect.(fx, w=fc_x) .* frect.(fy', w=fc_y)

    field_shifted = fftshift(fielmatrix(field))
    fft!(field_shifted)
    field_shifted .*= H
    ifft!(field_shifted)
    fieldmatrix(field) .= ifftshift(field_shifted)
end

export CoherentCirclePupil, CoherentEllipsePupil, CoherentRectanglePupil