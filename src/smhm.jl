# Functions relating to stellar-mass-halo-mass relations, typically derived from abundance matching.

"""
    flum(Mh::Real, mdef::String="vir")

Analytic approximations to the fraction of halos with mass definition `mdef` of mass `Mh` (in solar masses) that are luminous, based on Fig. 3 and 12 of Dooley 2017b. If `mdef="late z"`, then it is most similar to the purple dashed or solid green models in Figure 12. This was Annika's code. `mdef="200c"` is most similar to the Barber (blue) model in Figure 3. mdef="vir" is most similar to the yellow model in Figure 3. mdef="350c" is most similar to the red model in Figure 3. If `mdef="const"`, then just return `1`; e.g. 100% luminous fraction at all halo masses.

# Examples
```jldoctest
julia> flum(1e10, "late z") ≈ 0.9999938558253978
true

julia> flum(1e10, "200c") ≈ 0.9994793959258891
true

julia> flum(1e10, "vir") ≈ 0.9878715650157257
true

julia> flum(1e10, "350c") ≈ 0.9933071490757152
true

julia> flum(1e10, "const") ≈ 1.0
true
```
"""
function flum(Mh::Real, mdef::String="vir")
    if mdef == "late z"
        return 0.5 * (tanh( 3*(log10(Mh) - 8)) + 1)
    elseif mdef == "200c"
        return 0.5 * (tanh( 2.7*(log10(Mh) - 8.6)) + 1)
    elseif mdef == "vir"
        return 0.5 * (tanh( 2*(log10(Mh) - 8.9)) + 1)
    elseif mdef == "350c"
        return 0.5 * (tanh( 2*(log10(Mh) - 8.75)) + 1)
    elseif mdef == "const"
        return one(Mh)
    else
        error(mdef*" is not a valid string for `mdef`")
    end
end


"""
    Moster2013(Mh, z, M10=11.59, M11=1.195, N10=0.0351, N11=-0.0247, B10=1.376, B11=-0.826, gamma10=0.608, gamma11=0.329)
    
Given a halo mass `Mh` (M200c in solar masses) and a redshift `z`, return the median stellar mass from the [Moster et al. 2013](http://adsabs.harvard.edu/abs/2013MNRAS.428.3121M) empirical model. If `Mh` is a `Unitful.Mass`, will return a `Unitful.Mass`. 

# Examples
```jldoctest
julia> Moster2013( 1e12, 0.0 ) ≈ 3.427515185606664e10
true

julia> Moster2013( 1e12 * UnitfulAstro.Msun, 0.0 ) ≈ 3.427515185606664e10 * UnitfulAstro.Msun
true
```

# Notes
Moster et al. 2013 utilized the Millennium simulation ([Springel et al. 2005](https://ui.adsabs.harvard.edu/abs/2005Natur.435..629S/abstract)) for high-mass halos and the Millennium-II simulation ([Boylan-Kolchin et al. 2009](https://ui.adsabs.harvard.edu/abs/2009MNRAS.398.1150B/abstract)) for low-mass halos. At low redshift, they use the stellar mass function of [Li & White 2009](https://ui.adsabs.harvard.edu/abs/2009MNRAS.398.2177L) derived from SDSS DR7, ranging from stellar masses of 10^8.3 to 10^11.7 solar masses. This is supplemented by the low-mass stellar mass function from [Baldry et al. 2008](https://ui.adsabs.harvard.edu/abs/2008MNRAS.388..945B) for stellar masses between 10^7.4 to 10^8.3 solar masses. The high-redshift galaxy sample covers redshifts 0 < z < 4, using data from [Perez-Gonzalez et al. 2008](https://ui.adsabs.harvard.edu/abs/2008ApJ...675..234P) for high-mass galaxies and [Santini et al. 2012](https://ui.adsabs.harvard.edu/abs/2012A&A...538A..33S) for low-mass galaxies. See Figure 6 in [Moster et al. 2013](http://adsabs.harvard.edu/abs/2013MNRAS.428.3121M) for a visualization of the adopted stellar mass functions. 
"""
function Moster2013(Mh, z, M10=11.59, M11=1.195, N10=0.0351, N11=-0.0247, B10=1.376, B11=-0.826, gamma10=0.608, gamma11=0.329)
    #     M10, M11, N10, N11, B10, B11, gamma10, gamma11 = 11.59, 1.195, 0.0351, -0.0247, 1.376, -0.826, 0.608, 0.329
    zs = (z/(z+1))
    N=N10+N11*zs
    B=B10+B11*zs
    gamma=gamma10+gamma11*zs
    M1=exp10(M10+M11*zs)
    return 2 * Mh * N * ((((Mh/M1)^-B) + (Mh/M1)^gamma)^-1)
end
Moster2013(Mh::u.Mass, z, args...) = Moster2013(u.ustrip(ua.Msun,Mh), z, args...) * ua.Msun


"""
    Behroozi2013(Mh, z)

Given a halo mass `Mh` (specifically the peak virial halo mass in solar masses, as defined in [Bryan & Norman 1998](http://adsabs.harvard.edu/abs/1998ApJ...495...80B)) and a redshift `z`, return the median stellar mass from the [Behroozi et al. 2013](https://ui.adsabs.harvard.edu/abs/2013ApJ...770...57B/abstract) empirical model. If `Mh` is a `Unitful.Mass`, will return a `Unitful.Mass`.

# Examples
```jldoctest
julia> Behroozi2013( 1e12, 0.0 ) ≈ 2.6798246456860065e10
true

julia> Behroozi2013( 1e12 * UnitfulAstro.Msun, 0.0 ) ≈ 2.6798246456860065e10 * UnitfulAstro.Msun
true
```
"""
function Behroozi2013(Mh, z)
    a = 1 / (1 + z) # scale factor
    ν = exp(-4*a^2)
    ϵ = exp10( -1.777 + (-0.006 * (a-1))*ν - 0.119*(a-1) )
    M1 = exp10( 11.514 + (-1.793 * (a-1) + (-0.251*z))*ν )
    α = -1.412 + (0.731*(a-1))*ν
    δ = 3.508 + (2.608*(a-1) + (-0.043*z))*ν
    γ = 0.316 + (1.319 * (a-1) + 0.279*z)*ν
    MhICL = 12.515 + (-2.503 * (a-1))
    ρhalf = 0.799
    f(x) = -log10(exp10(α*x) + 1) + δ*(log10(1+exp(x)))^γ / (1+exp(exp10(-x)))
    return exp10( log10(ϵ*M1) + f(log10(Mh/M1)) - f(0.0))
end
Behroozi2013(Mh::u.Mass, z, args...) = Behroozi2013(u.ustrip(ua.Msun,Mh), z, args...) * ua.Msun


""" 
    GK14(Mh, z)

Given a halo mass `Mh` (specifically the peak virial halo mass in solar masses, as defined in [Bryan & Norman 1998](http://adsabs.harvard.edu/abs/1998ApJ...495...80B)) and a redshift `z`, return the median stellar mass from the [Garrison-Kimmel et al. 2014](https://ui.adsabs.harvard.edu/abs/2014MNRAS.438.2578G/abstract) empirical model. If `Mh` is a `Unitful.Mass`, will return a `Unitful.Mass`. This is the same model as [`Behroozi2013`](@ref), but with a shallower low-mass slope (α). See Figures 9 and 10. 

# Examples
```jldoctest
julia> GK14( 1e12, 0.0 ) ≈ 2.894144183405201e10
true

julia> GK14( 1e12 * UnitfulAstro.Msun, 0.0 ) ≈ 2.894144183405201e10 * UnitfulAstro.Msun
true
```
"""
function GK14(Mh, z)
    a = 1 / (1 + z) # scale factor
    ν = exp(-4*a^2)
    ϵ = exp10( -1.777 + (-0.006 * (a-1))*ν - 0.119*(a-1) )
    M1 = exp10( 11.514 + (-1.793 * (a-1) + (-0.251*z))*ν )
    α = -1.92 + (0.731*(a-1))*ν
    δ = 3.508 + (2.608*(a-1) + (-0.043*z))*ν
    γ = 0.316 + (1.319 * (a-1) + 0.279*z)*ν
    MhICL = 12.515 + (-2.503 * (a-1))
    ρhalf = 0.799
    f(x) = -log10(exp10(α*x) + 1) + δ*(log10(1+exp(x)))^γ / (1+exp(exp10(-x)))
    return exp10( log10(ϵ*M1) + f(log10(Mh/M1)) - f(0.0))
end
GK14(Mh::u.Mass, z, args...) = GK14(u.ustrip(ua.Msun,Mh), z, args...) * ua.Msun


""" 
    GK17_field(Mh, z, nu=-0.2)

Given a halo mass `Mh` (specifically the peak virial halo mass in solar masses, as defined in [Bryan & Norman 1998](http://adsabs.harvard.edu/abs/1998ApJ...495...80B)) and a redshift `z`, return the median stellar mass for field galaxies from the [Garrison-Kimmel et al. 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.464.3108G/abstract) empirical model. If `Mh` is a `Unitful.Mass`, will return a `Unitful.Mass`. This uses the model of [`Behroozi2013`](@ref) at high masses, but has a variable low-mass slope (α) and scatter (σν). The argument `nu` is the slope of the scatter in the SMHM relation (ν in their Equation 3). The default is -0.2 as in [Dooley et al. 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.471.4894D/abstract) "An Observers Guide ...".

# Examples
```jldoctest
julia> GK17_field( 1e12, 0.0, -0.2 ) ≈ 2.9300795754861504e10
true

julia> GK17_field( 1e12 * UnitfulAstro.Msun, 0.0, -0.2 ) ≈ 2.9300795754861504e10 * UnitfulAstro.Msun
true
```
"""
function GK17_field(Mh, z, nu=-0.2)
    a = 1 / (1 + z) # scale factor
    ν = exp(-4*a^2)
    ϵ = exp10( -1.777 + (-0.006 * (a-1))*ν - 0.119*(a-1) )
    M1 = exp10( 11.514 + (-1.793 * (a-1) + (-0.251*z))*ν )
    Mh>M1 ? (σν = 0.2) : (σν = 0.2 + nu * (log10(Mh) - log10(M1)))
    # σν = ifelse(Mh>M1, 0.2, 0.2 + nu * (log10(Mh) - log10(M1)))
    α = -(0.24 * σν^2 + 0.16*σν + 1.99) + (0.731*(a-1))*ν
    δ = 3.508 + (2.608*(a-1) + (-0.043*z))*ν
    γ = 0.316 + (1.319 * (a-1) + 0.279*z)*ν
    MhICL = 12.515 + (-2.503 * (a-1))
    ρhalf = 0.799
    # for some reason, need to provide α here to avoid allocations
    f(x) = -log10(exp10(α*x) + 1) + δ*(log10(1+exp(x)))^γ / (1+exp(exp10(-x)))
    return exp10( log10(ϵ*M1) + f(log10(Mh/M1)) - f(0.0))
end
GK17_field(Mh::u.Mass, z, args...) = GK17_field(u.ustrip(ua.Msun,Mh), z, args...) * ua.Msun


""" 
    GK17_sat(Mh, z, nu=-0.2)

Given a halo mass `Mh` (specifically the peak virial halo mass in solar masses, as defined in [Bryan & Norman 1998](http://adsabs.harvard.edu/abs/1998ApJ...495...80B)) and a redshift `z`, return the median stellar mass for satellite galaxies from the [Garrison-Kimmel et al. 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.464.3108G/abstract) empirical model. If `Mh` is a `Unitful.Mass`, will return a `Unitful.Mass`. This uses the model of [`Behroozi2013`](@ref) at high masses, but has a variable low-mass slope (α) and scatter (σν). The argument `nu` is the slope of the scatter in the SMHM relation (ν in their Equation 3). The default is -0.2 as in [Dooley et al. 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.471.4894D/abstract) "An Observers Guide ...".

# Examples
```jldoctest
julia> GK17_sat( 1e12, 0.0, -0.2 ) ≈ 2.860110341918958e10
true

julia> GK17_sat( 1e12 * UnitfulAstro.Msun, 0.0, -0.2 ) ≈ 2.860110341918958e10 * UnitfulAstro.Msun
true
```
"""
function GK17_sat(Mh, z, nu=-0.2)
    a = 1 / (1 + z) # scale factor
    ν = exp(-4*a^2)
    ϵ = exp10( -1.777 + (-0.006 * (a-1))*ν - 0.119*(a-1) )
    M1 = exp10( 11.514 + (-1.793 * (a-1) + (-0.251*z))*ν )
    Mh>M1 ? (σν = 0.2) : (σν = 0.2 + nu * (log10(Mh) - log10(M1)))
    # σν = ifelse(Mh>M1, 0.2, 0.2 + nu * (log10(Mh) - log10(M1)))
    α = -(0.14 * σν^2 + 0.14*σν + 1.79) + (0.731*(a-1))*ν
    δ = 3.508 + (2.608*(a-1) + (-0.043*z))*ν
    γ = 0.316 + (1.319 * (a-1) + 0.279*z)*ν
    MhICL = 12.515 + (-2.503 * (a-1))
    ρhalf = 0.799
    f(x) = -log10(exp10(α*x) + 1) + δ*(log10(1+exp(x)))^γ / (1+exp(exp10(-x)))
    return exp10( log10(ϵ*M1) + f(log10(Mh/M1)) - f(0.0))
end
GK17_sat(Mh::u.Mass, z, args...) = GK17_sat(u.ustrip(ua.Msun,Mh), z, args...) * ua.Msun
