module FastOptics

import FFTW: fft!, fft, ifft!, ifft,  fftshift, ifftshift, ifftshift!, fftfreq, Frequencies
import ImageTransformations: imresize
import Base: size, step, convert, range

# abstract types
include("abstracttypes.jl")

# types
include("types.jl")

# fields
include("fields/field.jl")

# utilities
include("utils/utils.jl")

# functions
include("functions/functions.jl")

# apertures
include("elements/elements.jl")

# imaging
include("imaging/imaging.jl")

# propagation methods
include("propagations/propagation.jl")

end
