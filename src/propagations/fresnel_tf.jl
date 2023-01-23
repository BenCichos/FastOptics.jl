struct FresnelTF <: PropagationFunction end

(::Type{FresnelTF})(field::LightField, distance::T, with_steps::Bool=false) where {T <: AbstractFloat} = fresnel_tf(field, distance, with_steps)

function fresnel_tf(field::LightField, distance::T, with_steps::Bool=false) where {T <: AbstractFloat} 
    field_space = fourier_space(field)
    lambda = wavelength(field)
    
    with_steps ? stepwise_fresnel_tf(field, field_space, lambda, distance) : fresnel_tf(field, field_space, lambda, distance)
end

function stepwise_fresnel_tf(field::LightField, field_space::Matrix{T}, lambda::T, distance::T) where {T <: AbstractFloat}
    steps, step_size = fresnel_tf_steps(field, distance)
    H = fresnel_transfer_function(field_space, lambda, step_size)
    
    field_shifted = fftshift(matrix(field))
    fft!(field_shifted)

    for _ in 1:steps
        field_shifted .*= H
    end    

    ifft!(field_shifted)
    ifftshift!(matrix(field), field_shifted)
end

function fresnel_tf(field::LightField, fourier_space::Matrix{T}, lambda::T, distance::T) where {T <: AbstractFloat}
    field_shifted = fftshift(matrix(field))
    fft!(field_shifted)
    field_shifted .*= fresnel_transfer_function(fourier_space, lambda, distance)
    ifft!(field_shifted)
    ifftshift!(matrix(field), field_shifted)
end

function fresnel_transfer_function(fourier_space::Matrix{T}, lambda::T, distance::T) where {T <: AbstractFloat}
    exp.((im*pi*lambda*distance) * fourier_space)
end

function fresnel_tf_steps(field::LightField, distance::T) where {T <: AbstractFloat} 
    x_critical, y_critical = critical_sampling(field)
    
    critical = max(x_critical, y_critical)
    steps = ceil(Int, distance / critical)
    step_size = distance / steps

    steps, step_size
end

export FresnelTF