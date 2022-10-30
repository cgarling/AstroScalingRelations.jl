module AstroScalingRelations

import Unitful
const u = Unitful
import UnitfulAstro
const ua = UnitfulAstro
import LambertW: lambertw

include("smhm.jl")
include("outflows.jl")
include("gas_mass.jl")

export flum, Moster2013, Behroozi2013, GK14, GK17_field, GK17_sat  # Exports from smhm.jl
export Muratov2015, Christensen2016 # Exports from outflows.jl
export Papastergis2012, Bradford2015, Scoville2017, Wang2016_DHI, Wang2016_rho_rs # Exports from gas_mass.jl

end # module
