##################
# Abstract Types #
##################


# Light Field
abstract type LightField end
abstract type MonochromaticLightField <: LightField end
abstract type PolychromaticLightField <: LightField end

export LightField, MonochromaticLightField, PolychromaticLightField


# Propagation Function
abstract type PropagationFunction end

export PropagationFunction


# Optical Element
abstract type AbstractOpticalElement end

# Aperture
abstract type AbstractAperture <: AbstractOpticalElement end
abstract type FunctionalAperture <: AbstractAperture end
abstract type GratingAperture <: FunctionalAperture end

export Aperture, FunctionAperture, GratingAperture
# Pupil

abstract type AbstractPupil <: AbstractOpticalElement end
abstract type AbstractCoherentPupil <: AbstractPupil end 
abstract type AbstractIncoherentPupil <: AbstractPupil end

export AbstractPupil, AbstractCoherentPupil, AbstractIncoherentPupil

