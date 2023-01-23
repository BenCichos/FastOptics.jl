struct FresnelTwoStep <: PropagationFunction end

(::Type{FresnelTwoStep})(field::LightField, distance::T, scale::T) where {T <: AbstractFloat} = fresnel_two_step(field, distance, scale)

function fresnel_two_step(field::LightField, distance::T, scale::T) where {T <: AbstractFloat}
    x_size, y_size = size(field)
    x_space, y_space = space(field)
    x_step, y_step = step(field)
    x_real, y_real = realsize(field)
    x_real_scaled, y_real_scaled = realsize(field) .* scale
    lambda = wavelength(field)
    k = wavenumber(field)
    
    fx = fftfreq(x_size, x_size / x_real)
    fy = fftfreq(y_size, y_size / y_real)
    
    matrix(field) .*= exp.(im*k/(2*distance*x_real)*(x_real-x_real_scaled)*(x_space.^2 .+ y_space'.^2))

    field_shifted = fftshift(matrix(field))    
    fft!(field_shifted)
    field_shifted .*= exp.(-im*pi*lambda*distance*x_real/x_real_scaled * (fx.^2 .+ fy'.^2))
    ifft!(field_shifted)
    matrix(field) .= ifftshift(field_shifted)

    x_step_scaled, y_step_scaled = (x_real_scaled, y_real_scaled) ./ (x_size, y_size)
    x_space_scaled = range(start=-x_real_scaled/2, step=x_step_scaled, length=x_size) 
    y_space_scaled = range(start=-y_real_scaled/2, step=y_step_scaled, length=y_size)
    realsize!(field, (x_real_scaled, y_real_scaled))
    space!(field, (x_space_scaled, y_space_scaled))

    normalise = (x_real_scaled * x_step * y_step) / (x_real * x_step_scaled * y_step_scaled)
    matrix(field) .*= normalise * fresnel_two_step_transfer_function(x_space_scaled, y_space_scaled, k, distance, x_real, x_real_scaled) 
end

function fresnel_two_step_transfer_function(x_space::AbstractVector, y_space::AbstractVector, k::AbstractFloat, distance::AbstractFloat, x_real::AbstractFloat, x_real_scaled::AbstractFloat)
    exp.(-im*k/(2*distance*x_real_scaled) * (x_real - x_real_scaled) * (x_space.^2 .+ y_space'.^2))
end

export FresnelTwoStep