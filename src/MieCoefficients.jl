"Computes Mie or Generalized Mie coefficients for a given object"
module MieCoefficients

# Dependencies
using MyPhysicalConstants, OpticalProperties, MultiLayerNFRHT, SpecialFunctions, MyPhysicalConstants

# Import for overloading functions
import OpticalProperties: permittivity, refractive_index

# export functions
export mie_coefficients, efficiencies, size_parameter, radius

# export types
export AbstractObject, Sphere, Cylinder


include("spheres.jl")
include("spherical_cavities.jl")
include("cylinders.jl")


end