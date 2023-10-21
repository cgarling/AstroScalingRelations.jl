module AstroScalingRelations

import Unitful
const u = Unitful
import UnitfulAstro
const ua = UnitfulAstro
import LambertW: lambertw

include("smhm.jl")
include("outflows.jl")
include("gas_mass.jl")

####################################
# Galaxy sizes

"""
    galaxy_size(Rvir; Aᵣ=0.02, fₚfₖ=0.78) = Aᵣ * fₚfₖ * Rvir
    galaxy_size(Mh, ρthresh; Aᵣ=0.02, fₚfₖ=0.78) = Aᵣ * fₚfₖ * cbrt( 3 * Mh / (4π * ρthresh) )

Scaling relation for galaxy half-light radius with the virial radius of the dark matter halo such that

```math
\\begin{aligned}
ρ_\\text{thresh} &= ρ_c(z) \\, Δ_\\text{vir} \\newline
R_\\text{vir} &= \\sqrt[3]{ \\frac{3 \\, \\text{M}_h}{4π \\; ρ_c(z)  \\; Δ_\\text{vir}}} \\newline
R_h(\\text{M}_h, z) &= A_r \\; f_p f_k \\; R_\\text{vir}
\\end{aligned}
```

The units of the returned ``R_h`` are the same as the provided ``R_\\text{vir}``, or if using the `(Mh, ρthresh)` signature, the unit of `cbrt( 3 * Mh / (4π * ρthresh) )`.

Default keyword arguments ``Aᵣ=0.02`` from [Jiang et al. 2019](https://doi.org/10.1093/mnras/stz1952) and ``fₚfₖ=0.78`` for the projection correction for spheroidal systems from [Somerville et al. 2018](https://doi.org/10.1093/mnras/stx2040).

# Examples
```jldoctest
julia> galaxy_size(1e12, 12000; Aᵣ=0.02, fₚfₖ=0.78) ≈ 4.227023347086455
true
```
"""
galaxy_size(Rvir; Aᵣ=0.02, fₚfₖ=0.78) = Aᵣ * fₚfₖ * Rvir
galaxy_size(Mh, ρthresh; Aᵣ=0.02, fₚfₖ=0.78) = Aᵣ * fₚfₖ * cbrt( 3 * Mh / (4π * ρthresh) )

####################################
# Exports

export flum, Moster2013, Behroozi2013, GK14, GK17_field, GK17_sat  # Exports from smhm.jl
export Muratov2015, Christensen2016 # Exports from outflows.jl
export Papastergis2012, Bradford2015, Scoville2017, Wang2016_DHI, Wang2016_rho_rs # Exports from gas_mass.jl
export galaxy_size # Exports from AstroScalingRelations.jl main file

end # module
