# Functions relating to gas outflows from galaxies due to stellar feedback; particularly mass-loading factors.

"""
    Muratov2015(Mstar)
    Muratov2015(Vvir, z)

The model for the mass-loading factor η from [Muratov et al. 2015](http://adsabs.harvard.edu/abs/2015MNRAS.454.2691M), based on FIRE-1 simulations. Two parameterizations are allowed; one which is based on the galaxy's stellar mass (`Mstar`) in solar masses and one which is based on the circular velocity at the virial radius (`Vvir`) of the galaxy's dark matter halo in km/s and the redshift of evaluation (`z`). Both parameterizations can be `Unitful` quantities. See also [`Pandya2021`](@ref). 

# Examples
```jldoctest
julia> Muratov2015( 1e10 ) ≈ 3.6
true

julia> Muratov2015( 1e10 * UnitfulAstro.Msun ) ≈ 3.6
true

julia> Muratov2015( 40, 1 ) ≈ 26.135392171726632
true

julia> Muratov2015( 40 * Unitful.km / Unitful.s, 1 ) ≈ 26.135392171726632
true
```
"""
Muratov2015(Mstar) = 3.6 * (Mstar/1e10)^(-0.35)
Muratov2015(Mstar::u.Mass) = Muratov2015(u.ustrip(ua.Msun, Mstar))
function Muratov2015(Vvir, z)
    α = ifelse(Vvir<60, -3.2, -1)
    return 2.9 * (1 + z)^1.3 * (Vvir/60)^α
end
Muratov2015(Vvir::u.Velocity, z) = Muratov2015(u.ustrip(u.km/u.s, Vvir), z)

"""
    Pandya2021(Mstar)
    Pandya2021(Vvir, z)

The model for the mass-loading factor η from [Pandya et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.508.2979P/abstract), based on FIRE-2 simulations. Two parameterizations are allowed; one which is based on the galaxy's stellar mass (`Mstar`) in solar masses and one which is based on the circular velocity at the virial radius (`Vvir`) of the galaxy's dark matter halo in km/s and the redshift of evaluation (`z`). Originally calculated for redshift ``0 ≤ z ≤ 4`` but will extrapolate outside this range. Both parameterizations can be `Unitful` quantities.

# Examples
```jldoctest
julia> isapprox( Pandya2021( 1e10 ), 0.6309573444801929; rtol=1e-7 )
true

julia> isapprox( Pandya2021( 1e10 * UnitfulAstro.Msun ), 0.6309573444801929; rtol=1e-7 )
true

julia> all( isapprox.( Pandya2021.( 40, (0.0, 1.0, 3.0) ), (12.977828436995544, 21.609465380663718, 24.88169815959356); rtol=1e-7 ) )
true

julia> isapprox( Pandya2021( 40 * Unitful.km / Unitful.s, 1 ), 21.609465380663718; rtol=1e-7 )
true
```
"""
Pandya2021(Mstar) = exp10(4.3) * Mstar^-0.45
Pandya2021(Mstar::u.Mass) = Pandya2021(u.ustrip(ua.Msun, Mstar))
function Pandya2021(Vvir, z)
    @assert z ≥ 0
    if z ≤ 0.5
        return exp10(6.4) * Vvir^-3.3
    elseif z ≤ 2
        return exp10(5.5) * Vvir^-2.6
    else
        return exp10(4.6) * Vvir^-2
    end
end
Pandya2021(Vvir::u.Velocity, z) = Pandya2021(u.ustrip(u.km/u.s, Vvir), z)

"""
    Christensen2016(Vcirc)

The model for the mass-loading factor η from [Christensen et al. 2016](http://adsabs.harvard.edu/abs/2016ApJ...824...57C), based on simulations run with the GASOLINE hydrodynamical model. The model requires the maximum circular velocity (`Vcirc`) of the galaxy's dark matter halo in km/s. `Vcirc` can be a `Unitful.Velocity`. Note that the normalization of the relation is not given explicitly in the paper, only the scaling is; the prefactor here was estimated from their Figure 11. 

# Examples
```jldoctest
julia> Christensen2016( 40 ) ≈ 3.0028111668298063
true

julia> Christensen2016( 40 * Unitful.km / Unitful.s ) ≈ 3.0028111668298063
true
```
"""
Christensen2016(Vcirc) = 10047.546 * Vcirc^-2.2
Christensen2016(Vcirc::u.Velocity) = Christensen2016(u.ustrip(u.km/u.s, Vcirc))
