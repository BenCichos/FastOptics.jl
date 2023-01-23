struct ImageAperture{T, S} <: AbstractAperture where {T <: AbstractFloat, S <: AbstractFloat}
    image::Matrix{T}
    size::NTuple{2, Union{S, Missing}}

    function ImageAperture(image::Matrix{T}, size::NTuple{2, Union{S, Missing}}) where {T <: AbstractFloat, S <: AbstractFloat}
        new{T, S}(image, size)
    end

end

function ImageAperture(image::Matrix{T}) where {T <: AbstractFloat}
    ImageAperture(image, (missing, missing))
end

apertureimage(ia::ImageAperture) = ia.image
imagesize(ia::ImageAperture) = size(ia.image)
size(ia::ImageAperture) = ia.size
width(ia::ImageAperture) = ia.size[1]
height(ia::ImageAperture) = ia.size[2]

export ImageAperture
export apertureimage, imagesize, size, width, height