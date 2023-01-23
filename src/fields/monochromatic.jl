#######################
# Monochromatic Field #
#######################

mutable struct MonochromaticField{T, S, R, F} <: MonochromaticLightField where {T <: AbstractFloat, S <: AbstractFloat, R <: AbstractRange, F <: Frequencies}
    matrix::Matrix{Complex{T}}
    wavelength::T
    wavenumber::T
    realsize::NTuple{2, S}
    range::NTuple{2, R}
    space::Matrix{T}
    fourier_range::NTuple{2, F}
    fourier_space::Matrix{T}

    function MonochromaticField(matrix::Matrix{Complex{T}}, wavelength::T, wavenumber::T, realsize::NTuple{2, S}, range::NTuple{2, R}, space::Matrix{T}, fourier_range::NTuple{2, F}, fourier_space::Matrix{T}) where {T <: AbstractFloat, S <: AbstractFloat, R <: AbstractRange, F <: Frequencies }
        new{T, S, R, F}(matrix, wavelength, wavenumber, realsize, range, space, fourier_range, fourier_space)
    end
end

function MonochromaticField(matrix::Matrix{Complex{T}}, wavelength::T, realsize::Tuple{S, S}) where {T <: AbstractFloat, S <: AbstractFloat}
    wavenumber = 2 * pi / wavelength

    x_size, y_size = size(matrix)
    x_real, y_real = realsize
    x_step, y_step = x_real / x_size, y_real / y_size
    x_range = range(start=-x_real/2, step=x_step, length=x_size)
    y_range = range(start=-y_real/2, step=y_step, length=y_size)
    
    fx_range = fftfreq(x_size, x_size / x_real)
    fy_range = fftfreq(y_size, y_size / y_real)

    space = x_range.^2 .+ y_range'.^2
    fourier_space = fx_range.^2 .+ fy_range'.^2
    
    MonochromaticField(matrix, wavelength, wavenumber, realsize, (x_range, y_range), space, (fx_range, fy_range), fourier_space)
end

function MonochromaticField(wavelength::T, range::NTuple{2, R}) where {T <: AbstractFloat, R <: StepRangeLen}
    wavenumber = 2 * pi / wavelength
    
    x_range, y_range = range
    x_size, y_size = length(x_range), length(y_range)
    x_real, y_real = x_range[end] - x_range[1], y_range[end] - y_range[1]
    
    fx_range = fftfreq(x_size, x_size / x_real)
    fy_range = fftfreq(y_size, y_size / y_real)
    
    space = x_range .+ y_range'
    fourier_space = fx_range .+ fy_range' 
    
    matrix = ones(ComplexF64, (x_size, y_size))

    MonochromaticField(matrix, wavelength, wavenumber, realsize, range, space, (fx_range, fy_range), fourier_space)
end

matrix(f::MonochromaticField) = f.matrix
size(f::MonochromaticField) = size(f.matrix)
wavelength(f::MonochromaticField) = f.wavelength
wavenumber(f::MonochromaticField) = f.wavenumber
realsize(f::MonochromaticField) = f.realsize
range(f::MonochromaticField) = f.range
space(f::MonochromaticField) = f.space
fourier_range(f::MonochromaticField) = f.fourier_range
fourier_space(f::MonochromaticField) = f.fourier_space
step(f::MonochromaticField) = step.(range(f))
fourier_step(f::MonochromaticField) = step.(fourier_range(f))

matrix!(f::MonochromaticField, new::Matrix{ComplexF64}) = setfield!(f, :matrix, new )
realsize!(f::MonochromaticField, new::NTuple{2, T}) where {T <: AbstractFloat} = setfield!(f, :realsize, new)
space!(f::MonochromaticField, new::Matrix{T}) where {T <: AbstractFloat} = setfield!(f, :space, new)

export MonochromaticField
export matrix, size, wavelength, wavenumber, realsize, range, space, fourier_range, fourier_space, step, fourier_step
export matrix!, realsize!, space!
