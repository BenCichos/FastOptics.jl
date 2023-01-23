function _crop_image_one!(field::LightField, aperture_image::Matrix, size::Int, aperture_size::Int; crop_dimension::CropDimension)
    diff = size - aperture_size
    start_offset = end_offset = floor(Int, abs(diff) / 2)
    iseven(diff) && (start_offset += 1)

    if diff > 0
        _cropped_view_one(matrix(field), crop_dimension, start_offset, size - end_offset - 1) .*= aperture_image
        _cropped_view_one(matrix(field), crop_dimension, 1, start_offset - 1) .= 0.0
        _cropped_view_one(matrix(field), crop_dimension, size - end_offset, size) .= 0.0
    else
        _cropped_view_one(aperture_image, crop_dimension, start_offset, aperture_size-end_offset)
    end
end

function _crop_image_two!(field::LightField, aperture_image::Matrix, size::Tuple{Int, Int}, aperture_size::Tuple{Int, Int})
    diff = size .- aperture_size
    start_offset = end_offset = floor.(Int, abs.(diff) ./ 2)
    is_even = iseven.(diff)
    start_offset = (start_offset[1] + (is_even[1] ? 0 : 1), start_offset[2] + (is_even[2] ? 0 : 1))

    if diff[1] > 0 && diff[2] > 0 
        _cropped_view_two(matrix(field), start_offset, size .- end_offset .- 1) .*= aperture_image
        _cropped_view_two(matrix(field), (1, 1), start_offset .- 1) .= 0.0
        _cropped_view_two(matrix(field), size .- end_offset, size) .= 0.0
    elseif diff[1] > 0 && diff[2] < 0
        _cropped_view_one(matrix(field), _crop_x, start_offset[1],  size[1] - end_offset[1] - 1) .*= _cropped_view_one(aperture_image, _crop_y, start_offset[2], aperture_size[2] - end_offset[2] - 1)
        _cropped_view_one(matrix(field), _crop_x, 1, start_offset[1] - 1) .= 0.0
        _cropped_view_one(matrix(field), _crop_x, size[1] - end_offset[1],  size[1]) .= 0.0
    elseif diff[1] < 0 && diff[2] > 0
        _cropped_view_one(matrix(field), _crop_y, start_offset[2],  size[2] - end_offset[2] - 1) .*= _cropped_view_one(aperture_image, _crop_x, start_offset[1], aperture_size[1] - end_offset[1] - 1)
        _cropped_view_one(matrix(field), _crop_y, 1, start_offset[2] - 1) .= 0.0
        _cropped_view_one(matrix(field), _crop_y, size[2] - end_offset[2],  size[2]) .= 0.0
    else
        matrix(field) .*= _cropped_view_two(aperture_image, start_offset, aperture_size .- end_offset .- 1)
    end
end

_cropped_view_one(matrix::Matrix, cd::CropDimension, first::Int, last::Int) = (cd === _crop_x) ? view(matrix, first:last, :) : view(matrix, :, first:last)
_cropped_view_two(matrix::Matrix, first::Tuple{Int, Int}, last::Tuple{Int, Int}) = view(matrix, first[1]:last[1], first[2]:last[2])