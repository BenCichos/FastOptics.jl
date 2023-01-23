####################
# Grating Aperture #
####################


# Methods For Subtypes of GratingAperture

function slitwidth end
function spacing end
function width end

spacing(ga::GratingAperture) = ga.spacing
width(ga::GratingAperture) = ga.width


# Cosine Wave Grating

struct CosineWaveGrating{T} <: GratingAperture
    fn::Function
    width::T
    spacing::T
    
    function CosineWaveGrating(width::T, spacing::T) where {T <: AbstractFloat}
        (width > 0) || throw(ArgumentError("The width of the grating has to be a positive number."))
        (spacing > 0) || throw(ArgumentError("The spacing of the grating has to be a positive number."))

        fn(field::LightField) = cosine_wave_grating(field=field, width=width, spacing=spacing)
        new{T}(fn, width, spacing)
    end
end

CosineWaveGrating(;width::T, spacing::T) where {T <: AbstractFloat} = CosineWaveGrating(width, spacing)

function cosine_wave_grating end

function cosine_wave_grating(field::LightField, width::Float64, spacing::Float64)
    x_range , y_range = range(field)
    period = 1 .- cos.(2 * pi * x_range / spacing)
    grating_one_dim = 1/2 * period .* frect.(x_range, w=width)
    grating_one_dim .* frect.(y_range', w=width)
end

cosine_wave_grating(;field::LightField, width::Float64, spacing::Float64) = cosine_wave_grating(field, width, spacing) 


# Methods on CosineWaveGrating

slitwidth(cwg::CosineWaveGrating) = cwg.spacing / 2


# SquareWaveGrating

struct SquareWaveGrating{T} <: FunctionalAperture
    fn::Function
    width::T
    spacing::T
    slitwidth::T
    
    function SquareWaveGrating(width::T, spacing::T, slitwidth::T) where {T <: AbstractFloat}
        (width > 0) || throw(ArgumentError("The width of the grating has to be positive number."))
        (spacing > 0) || throw(ArgumentError("The spacing of the grating has to be a positive number."))
        (slitwidth > 0 && slitwidth < spacing) || throw(ArgumentError("The slitwidth of the grating has to be less than the slit spacing"))

        fn(field::LightField) = square_wave_grating(field=field, width=width, spacing=spacing, slitwidth=slitwidth)
        new{T}(fn, width, spacing, slitwidth)
    end
end

SquareWaveGrating(;width::Float64, spacing::Float64, slitwidth::Float64) = SquareWaveGrating(width, spacing, slitwidth)

function square_wave_grating end

function square_wave_grating(;field::LightField, width::Float64, spacing::Float64, slitwidth::Float64)
    @assert width > 0 "The width of the grating has to be positive number"
    @assert spacing > 0 "The spacing of the grating has to be a positive number"
    @assert slitwidth > 0 "The slit width of the grating has to be a positive number"
    @assert spacing > slitwidth "The grating spacing has to be larger than the slit width"
    
    x_size, _ = size(field)
    x_real, _ = realsize(field)
    _, y_range = range(field)
    x_step, _ = step(field)
    
    # calculate grating spacing in terms of array indices
    number_of_slits = round(Int, width / spacing )
    spacing_array_space = round(Int, (spacing / x_real) * x_size )
    
    # calculate the slit width in terms of array indices
    slitwidth_array_space = round(Int, (slitwidth / x_real) * x_size)

    # check that the matrix has enough resolution to display grating
    @assert spacing_array_space > 0 "The spacing of the grating is too small. Either increase the size of the matrix or reduce the real size of the field."
    @assert slitwidth_array_space > 0 "The slit width of the grating is too small. Either increase the size of the matrix or reduce the real size of the field."
    
    # compute starting and end index of grating
    start_index = round(Int, (x_real - width) / (2 * x_step)) + 1
    end_index = start_index + (spacing_array_space * number_of_slits)
        
    # check that indices are within the bounds of the array
    (start_index < 1) && (start_index = 1)
    (end_index > x_size) && (end_index = x_size)

    # initialise bitvector to store one dimensional grating
    grating_line = falses(x_size)

    # create square wave grating
    @inbounds @simd for index = start_index:spacing_array_space:end_index-1
        grating_line[index: index+slitwidth_array_space] .= true
    end
    
    grating_line .* frect.(y_range', w=width)
end


# Methods for SquareWaveGrating

slitwidth(swg::SquareWaveGrating) = swg.slitwidth

export width, spacing, slitwidth
export CosineWaveGrating, cosine_wave_grating
export SquareWaveGrating, square_wave_grating