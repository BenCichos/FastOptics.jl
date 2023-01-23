struct FresnelIR <: PropagationFunction end

(::Type{FresnelIR})(field::LightField, distance::T) where {T <: AbstractFloat} = fresnel_ir(field, distance)

function fresnel_ir(field::LightField, distance::T) where {T <: AbstractFloat} 
    field_space = space(field)
    lambda = wavelength(field)
    k = wavenumber(field)

    space_shifted = fftshift(field_space)
    field_shifted = fftshift(matrix(field))
    
    fft!(field_shifted)
    field_shifted .*=  prod(step(field)) / (im*lambda*distance)
    field_shifted .*= fft!(fresnel_ir_function(space_shifted, distance, k))
    ifft!(field_shifted)

    ifftshift!(matrix(field), field_shifted)
end

@inline function fresnel_ir_function(space_shifted::Matrix{T}, distance::T, k::T) where {T <: AbstractFloat} 
    exp.( ( im*k / (2*distance) ) * space_shifted )
end

export FresnelIR