###############
# Propagation #
###############


# Fresnel Propagation Functions

## Fresnel Transfer Function - Nearfield Propagation
include("fresnel_tf.jl")

## Fresnel Impulse Response Function - Farfield Propagation
include("fresnel_ir.jl")

## Fresnel Two Step Transfer Function - Nearfield Propagation / Scaling
include("fresnel_twostep.jl")


# Fraunhofer Propagation Functions

include("fraunhofer.jl")

# Define propagate! methods

function propagate! end 

propagate!(field::LightField, distance::T) where {T <: AbstractFloat} = any(distance .> critical_sampling(field)) ? propagate!(field, distance, FresnelTF) : propagate!(field, distance, FresnelIR) 
propagate!(field::LightField, distance::T, scale::T) where {T <: AbstractFloat} =  FresnelTwoStep(field, distance, scale)
propagate!(field::LightField, distance::T, prop_fn::Type{<:PropagationFunction}) where {T <: AbstractFloat} = isa(prop_fn, FresnelTwoStep) ? prop_fn(field, distance, 1.0) : prop_fn(field, distance)


# Define propagate methods for each specific propagation functions

fresnel_tf_propagate!(field::LightField, distance::T, with_steps::Bool=false) where {T <: AbstractFloat} = FresnelTF(field, distance, with_steps)
fresnel_ir_propagate!(field::LightField, distance::T) where {T <: AbstractFloat} = propagate!(field, distance, FresnelIR)
fresnel_two_step_propagate!(field::LightField, distance::T, scale::T=1.0) where {T <: AbstractFloat} = propagate!(field, distance, scale)
fraunhofer_propagate!(field::LightField, distance::T) where {T <: AbstractFloat} = propagate!(field, distance, Fraunhofer)

export propagate!, fresnel_tf_propagate!, fresnel_ir_propagate!, frensel_two_step_propagate!, fraunhofer_propagate!


