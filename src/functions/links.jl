# Extend LinkedFunctions Package

fft!(l::LinkedFunction) = (l.is_fn1 = !l.is_fn1; l)
fft(l::LinkedFunction) = other(l)

export fft, fft!

# Link frect with fsinc
const linked_frect = @link frect fsinc
export linked_frect, @frect, @fsinc

@linkconvert @frect @fsinc function convert(x::Real; w::Real)
    (x,), (w = 1/w,)
end

@linkconvert @frect @fsinc function convert(x::Real; w::Real)
    (x,), (w = 1/w,)
end

const linked_frect2d = @link frect2d fsinc2d
const linked_fgaussian = @link fgaussian fgaussian
