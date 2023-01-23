###################
# OpticalElements #
###################

# Generic function for adding elements
function add! end 
export add!

#############
# Apertures #
#############

# Define generic functions that act as convenience methods to apply an aperture to a field

# Generic function that applies aperture to a given field
function aperture! end
export aperture!

# Genric function that applied aperture to a given field in fourier space
function fftaperture! end
export fftaperture!


########################
# Functional Apertures #
########################

include("rectangle.jl")
include("circular.jl")
include("grating.jl")

# Methods For Functional Apertures
fn(fa::FunctionalAperture) = fa.fn
(fa::FunctionalAperture)(args...; kwargs...) = fn(fa)(args...; kwargs...)

# Define add! method for apertures of type Functional Aperture
function add!(field::LightField, aperture::FunctionalAperture)
    @inline matrix(field) .*= fn(aperture)(field)
end

# Define aperture! method for apertures of type FunctionalAperture
function aperture!(field::LightField, aperture::FunctionalAperture)
    @inline add!(field, aperture)
end

# Define fftaperture! functions for FunctionalApertures

function fftaperture!(field::LightField, aperture::RectangleAperture)
    x_space, y_space = space(field)
    x_real, y_real = realsize(field)
    x_size, y_size = size(field)
    x_step, y_step = step(field)

    width_x, width_y = size(aperture)
    f_width_x = width_x / (x_real * x_step)
    f_width_y = width_y / (y_real * y_step)
    aperture = RectangleAperture(f_width_x, f_width_y)

    fx = fftshift(fftfreq(x_size, x_size / x_real))
    fy = fftshift(fftfreq(y_size, y_size / y_real))
    
    space!(field, (fx, fy))
    fft!(matrix(field))

    aperture!(field, aperture)

    ifft!(matrix(field))
    space!(field, (x_space, y_space))
end

function fftaperture!(field::LightField, aperture::CircleAperture)
    x_space, y_space = space(field)
    x_real, y_real = realsize(field)
    x_size, y_size = size(field)
    x_step, y_step = step(field)

    radius = radius(aperture)
    f_radius = radius / (x_real * x_step)
    aperture = CircleAperture(f_radius)

    fx = fftshift(fftfreq(x_size, x_size / x_real))
    fy = fftshift(fftfreq(y_size, y_size / y_real))
    
    space!(field, (fx, fy))
    fft!(matrix(field))

    aperture!(field, aperture)

    ifft!(matrix(field))
    space!(field, (x_space, y_space))
end

function fftaperture!(field::LightField, aperture::EllipseAperture)
    x_space, y_space = space(field)
    x_real, y_real = realsize(field)
    x_size, y_size = size(field)
    x_step, y_step = step(field)

    radius_x, radius_y = size(aperture)
    f_radius_x = radius_x / (x_real * x_step)
    f_radius_y = radius_y / (y_real * y_step)
    aperture = EllipseAperture(f_radius_x, f_radius_y)

    fx = fftshift(fftfreq(x_size, x_size / x_real))
    fy = fftshift(fftfreq(y_size, y_size / y_real))
    
    space!(field, (fx, fy))
    fft!(matrix(field))

    aperture!(field, aperture)

    ifft!(matrix(field))
    space!(field, (x_space, y_space))
end


##################
# Image Aperture #
##################


include("image.jl")

# Define add! method for apertures of type ImageAperture
function add!(field::LightField, aperture::ImageAperture)
    aperture_x_real, aperture_y_real = size(aperture)
    aperture_x_real_missing = ismissing(aperture_x_real)
    aperture_y_real_missing = ismissing(aperture_y_real)
    
    x_size, y_size = size(field)
    x_real, y_real = realsize(field)
    
    x_step, y_step = step(field)
    x_real_diff = isapprox(x_real, aperture_x_real, atol=x_step/2)
    y_real_diff = isapprox(y_real, aperture_y_real, atol=y_step/2)

    x_scaling = x_real / aperture_x_real
    y_scaling  = y_real / aperture_y_real
    
    if aperture_x_real_missing && aperture_y_real_missing
        aperture_image = apertureimage(aperture)
    elseif !aperture_x_real_missing && !aperture_y_real_missing
        if x_real_diff && y_real_diff
            aperture_image  = imresize(apertureimage(aperture), ratio=(x_scaling, y_scaling))
        elseif x_real_diff
            aperture_image = imresize(apertureimage(aperture), ratio=(x_scaling, 1))
        elseif y_real_diff
            aperture_image = imresize(apertureimage(aperture), ratio=(1, y_scaling))
        end
    elseif aperture_x_real_missing && y_real_diff
        aperture_image = imresize(apertureimage(aperture), ratio=(1, y_scaling))
    elseif aperture_y_real_missing && !x_real
        aperture_image = imresize(apertureimage(aperture), ratio=(x_scaling, 1))
    end

    aperture_x_size, aperture_y_size = size(aperture_image)

    x_equal = x_size == aperture_x_size
    y_equal = y_size == aperture_y_size

    if x_equal && y_equal
        # aperture is the same size as the field
        matrix(field) .*= apertureimage(aperture)
    elseif x_equal
        # x dimension of aperture is the same modify y dimension
        _crop_image_one!(field, aperture_image, y_size, aperture_y_size, crop_dimension=_crop_y)
    elseif y_equal
        # y dimension of aperture is the same modify x dimension
        _crop_image_one!(field, aperture_image, x_size, aperture_x_size, crop_dimension=_crop_x)
    else
        # both dimension of the aperture are different
       _crop_image_two!(field, aperture_image, size(field), size(aperture_image))
    end

    matrix(field)
end

# Define aperture! method for apertures of type ImageAperture
function aperture!(field::LightField, aperture::ImageAperture)
    @inline add!(field, aperture)
end


############
# ThinLens #
############

include("lens.jl")

# Define add! method for a thin lens
function add!(field::LightField, lens::ThinLens)
    fieldspace = space(field)
    k = wavenumber(field)
    matrix(field) .*= exp.(-im * k * (fieldspace) / (2 * focal_length(lens)))
end

# Define thinlens! method for a thin lens
function thinlens!(field::LightField, focal_length::T) where {T <:AbstractFloat}
    @inline add!(field, ThinLens(focal_length))
end


########
# Tilt #
########

include("tilt.jl")

# Define add! method for a tilt
function add!(field::LightField, tilt::Tilt)
    x_range, y_range = range(field)
    k = wavenumber(field)
    θ, α = angle(tilt)

    x_prime = x_range * cos(θ)
    y_prime = y_range * sin(θ)

    matrix(field) .*= exp.(im * k * (x_prime .+ y_prime') * tan(α))
end

# Define tilt! method for a tilt
function tilt!(field::LightField, theta::T, alpha::T) where {T <: AbstractFloat}
    @inline add!(field, Tilt(theta, alpha))
end
