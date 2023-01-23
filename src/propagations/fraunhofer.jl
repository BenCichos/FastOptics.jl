struct Fraunhofer <: PropagationFunction end

(::Type{Fraunhofer})(field::LightField, distance::T) where {T <: AbstractFloat} = fraunhofer(field, distance)

function fraunhofer(field::LightField, distance::T) where {T <: AbstractFloat} 
    field_space = space(field)
    lambda = wavelength(field)
    k = wavenumber(field)
    
    space_scaled = inv.(field_space) * lambda * distance
    real_scaled  = inv.(realsize(field)) .* (lambda * distance)

    space!(field, space_scaled)
    realsize!(field, real_scaled)
    
    field_shifted = fftshift(matrix(field))
    fft!(field_shifted)
    ifftshift!(matrix(field), field_shifted)
    matrix(field) .*= fraunhofer_function(space_scaled, lambda, distance, k)
end

function fraunhofer_function(space_scaled::Matrix{T},lambda::T, distance::T, k::T) where {T <: AbstractFloat}
    (1 / (im * lambda * distance)) * exp.((im*k / (2*distance)) * space_scaled )
end

export Fraunhofer