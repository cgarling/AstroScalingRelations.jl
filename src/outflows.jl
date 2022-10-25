# Functions relating to gas outflows from galaxies due to stellar feedback; particularly mass-loading factors.

"""
    Muratov2015(Mstar)
    Muratov2015(Vcirc, z)

The model for the mass-loading factor η from [Muratov et al. 2015](http://adsabs.harvard.edu/abs/2015MNRAS.454.2691M), based on FIRE simulations. Two parameterizations are allowed; one which is based on the galaxy's stellar mass (`Mstar`) in solar masses and one which is based on the maximum circular velocity (`Vcirc`) of the galaxy's dark matter halo in km/s and the redshift of evaluation (`z`). Both parameterizations can be `Unitful` quantities. See also Pandya et al. 2021.

# Examples
```jldoctest
julia> Muratov2015( 1e10 )
3.6

julia> Muratov2015( 1e10 * UnitfulAstro.Msun )
3.6

julia> Muratov2015( 40, 1 )
26.135392171726632

julia> Muratov2015( 40 * Unitful.km / Unitful.s, 1 )
26.135392171726632
```
"""
Muratov2015(Mstar) = 3.6 * (Mstar/1e10)^(-0.35)
Muratov2015(Mstar::u.Mass) = Muratov2015(u.ustrip(ua.Msun, Mstar))
function Muratov2015(Vcirc, z)
    α = ifelse(Vcirc<60, -3.2, -1)
    return 2.9 * (1 + z)^1.3 * (Vcirc/60)^α
end
Muratov2015(Vcirc::u.Velocity, z) = Muratov2015(u.ustrip(u.km/u.s, Vcirc), z)

"""
    Christensen2016(Vcirc)

The model for the mass-loading factor η from [Christensen et al. 2016](http://adsabs.harvard.edu/abs/2016ApJ...824...57C), based on simulations run with the GASOLINE hydrodynamical model. The model requires the maximum circular velocity (`Vcirc`) of the galaxy's dark matter halo in km/s. `Vcirc` can be a `Unitful.Velocity`. Note that the normalization of the relation is not given explicitly in the paper, only the scaling is; the prefactor here was estimated from their Figure 11. 

# Examples
```jldoctest
julia> Christensen2016( 40 )
3.0028111668298063

julia> Christensen2016( 40 * Unitful.km / Unitful.s )
3.0028111668298063
```
"""
Christensen2016(Vcirc) = 10047.546 * Vcirc^-2.2
Christensen2016(Vcirc::u.Velocity) = Christensen2016(u.ustrip(u.km/u.s, Vcirc))
