##########
# Fields #
##########

# Monochromatic Field
include("monochromatic.jl")

# Polychromatic Field

critical_sampling(field::LightField) = step(field) .* realsize(field) ./ wavelength(field)

export critical_sampling