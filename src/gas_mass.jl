# Functions relating to galactic gas masses and density profiles, mostly empirical.

"""
    Papastergis2012(Mstar)

Returns the mean neutral hydrogen gas mass of a galaxy with stellar mass `Mstar` in solar masses; this is the black line in Figure 19 of [Papastergis et al. 2012](http://adsabs.harvard.edu/abs/2012ApJ...759..138P). Based on a crossmatch between SDSS galaxies and sources detected in the blind ALFALFA α.40 survey. If `Mstar` is a `Unitful.Mass`, will return a `Unitful.Mass`.

# Examples
```jldoctest
julia> Papastergis2012(1e8) ≈ 2.0417379446695232e8
true

julia> Papastergis2012(1e8 * UnitfulAstro.Msun) ≈ 2.0417379446695232e8 * UnitfulAstro.Msun
true
```
"""
Papastergis2012(Mstar) = exp10(0.57 * log10(Mstar) + 3.75)
Papastergis2012(Mstar::u.Mass) = Papastergis2012(u.ustrip(ua.Msun,Mstar)) * ua.Msun


"""
    Bradford2015(Mstar)

Returns the median neutral hydrogen gas mass of an isolated galaxy with stellar mass `Mstar` in solar masses according to Equation in [Bradford et al. 2015](https://ui.adsabs.harvard.edu/abs/2015ApJ...809..146B/abstract). Uses coefficients from the first two rows of Table 3. These are derived for galaxies selected from SDSS DR8 with stellar masses between `exp10(7)` and `exp10(9.5)`.

# Examples
```jldoctest
julia> Bradford2015(1e8) ≈ 4.487453899331332e8
true

julia> Bradford2015(1e8 * UnitfulAstro.Msun) ≈ 4.487453899331332e8 * UnitfulAstro.Msun
true
```
"""
Bradford2015(Mstar) = Mstar < exp10(8.6) ? exp10(1.052*log10(Mstar) + 0.236) : exp10(0.461*log10(Mstar) + 5.3186) # 5.329) # 5.329 is the published value, but it is not continuous
Bradford2015(Mstar::u.Mass) = Bradford2015(u.ustrip(ua.Msun, Mstar)) * ua.Msun


"""
    Scoville2017(Mstar [, z, SFR])

Returns the median neutral hydrogen gas mass of a galaxy with stellar mass `Mstar` in solar masses at redshift `z` according to Equation 6 in [Scoville et al. 2017](http://adsabs.harvard.edu/abs/2017ApJ...837..150S). This does not include the specfic star formation rate dependence. These are derived from ALMA dust continuum emission (see Section 4) for well-studied galaxies in the COSMOS field between redshifts 0.3 and 4.5. If `Mstar` is a `Unitful.Mass`, will return a `Unitful.Mass`.

# Examples
```jldoctest
julia> Scoville2017(1e8) ≈ 1.7759037070772731e9
true

julia> Scoville2017(1e8 * UnitfulAstro.Msun) ≈  1.7759037070772731e9 * UnitfulAstro.Msun
true

julia> Scoville2017(1e8, 1.0) ≈ 6.357913365552342e9
true

julia> Scoville2017(1e8 * UnitfulAstro.Msun, 1.0) ≈ 6.357913365552342e9 * UnitfulAstro.Msun
true
```
"""
Scoville2017(Mstar) = 7.07e9 * (Mstar/1e10)^0.3
Scoville2017(Mstar::u.Mass) = Scoville2017(u.ustrip(ua.Msun,Mstar)) * ua.Msun
Scoville2017(Mstar, z) = 7.07e9 * (1 + z)^1.84 * (Mstar/1e10)^0.3
Scoville2017(Mstar::u.Mass, z) = Scoville2017(u.ustrip(ua.Msun,Mstar), z) * ua.Msun

###########################################################################################

"""
    Wang2016_DHI(MHI)

Returns the HI diameter, defined as the diameter at which the HI surface density equals 1 Msun/pc^2, in kiloparsecs given an HI mass follwing Equation 2 in [Wang et al. 2016](https://ui.adsabs.harvard.edu/abs/2016MNRAS.460.2143W/abstract). If `Mstar` is a `Unitful.Mass`, will return a `Unitful.Length`.

# Examples
```jldoctest
julia> Wang2016_DHI(1e8) ≈ 5.688529308438413
true

julia> Wang2016_DHI(1e8 * UnitfulAstro.Msun) ≈ 5.688529308438413 * UnitfulAstro.kpc
true
```
"""
Wang2016_DHI(MHI) = exp10(0.506 * log10(MHI) - 3.293)
Wang2016_DHI(MHI::u.Mass) = Wang2016_DHI(u.ustrip(ua.Msun,MHI)) * ua.kpc

"""
    (Σ₀, rs) = Wang2016_rho_rs(MHI)

Returns the central surface density `Σ₀` in Msun/pc^2 and the scale radius `rs` in parsecs that defines the surface density profile ``Σ(r) = Σ₀ \\exp (-r/rs)`` for a galaxy with an HI disk mass of `MHI` in solar masses, using the [Wang et al. 2016](https://ui.adsabs.harvard.edu/abs/2016MNRAS.460.2143W/abstract) relation between HI mass and HI disk diameter.

# Examples
```jldoctest
julia> all( Wang2016_rho_rs(1e6) .≈ (4.009796759878976, 199.2273166291845) )
true

julia> all( Wang2016_rho_rs(1e12) .≈ (7.045959041825033, 150293.43437081948) )
true
```
"""
function Wang2016_rho_rs(MHI)
    # Get the diameter at which the HI surface density equals 1 solmass / pc^2
    DHI = Wang2016_DHI(MHI) * 1e3 # Units of parsecs
    # We now have two constraints and two variables. 
    # Constraints are that the integral from 0 to ∞ of 2πr Σ0 exp(-r/rs) = MHI
    # This integral is analytic with result MHI = 2π Σ0 rs^2
    # and Σ0 exp(-(DHI/2)/rs) = 1 solmass / pc^2.
    # We'll solve first for Σ0, Σ0 = MHI / (2π rs^2) and substitute into the second
    # MHI / (2π rs^2) exp(-(DHI/2)/rs) = 1 solmass / pc^2
    # In Mathematica, we can run
    # Solve[MHI/(2 \[Pi]  rs^2) Exp[-(DHI/2)/rs] == 1, rs]
    # The resulting expression contains ProductLog, which is Mathematica's lambertw
    # These expressions can be further diffentiated with respect to MHI
    # so that they can be tracked as the gas mass changes. 
    # rs = -DHI/(4 * lambertw(-0.5 * DHI * sqrt(π/(2*MHI))))
    # C = 8*MHI/(DHI^2*ℯ^2*π)
    C = -0.5 * DHI * sqrt(π / (2 * MHI))
    if C > -1/ℯ
        # rs = -DHI/(4 * lambertw(-0.5 * DHI * sqrt(π/(2*MHI)),0))
        rs = -DHI / (4 * lambertw(C, 0))
    # elseif C == 1
    #     rs = DHI/4
    # else
    #     rs = -DHI/(4 * real(lambertw(-0.5 * DHI * sqrt(π/(2*MHI))+0im,-1)))
    else
        # rs = -DHI/(4 * (real(lambertwbp(C+1/ℯ+0im,0))-1))
        # rs = -DHI/(4 * lambertw(-1/ℯ+(-C-1/ℯ),-1))
        # rs = -DHI/(4 * (lambertwbp(-C-1/ℯ,-1)-1))
        ### rs=DHI/4 gives very similar answers to the alternate lambert branch but its not exact
        rs = DHI / 4
        # rs = DHI*0.26
    end
    # rs = -DHI/(4 * lambertw(-0.5 * DHI * sqrt(π/(2*MHI))))
    # rs = -DHI/(4 * real(lambertwbp(C+1/ℯ+0im,0)-1))
    Σ0 = MHI / (2π * rs^2)
    return (Σ0, rs) # return (rs, Σ0)
end
