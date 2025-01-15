""" `AbstractMZR{T <: Real}`: abstract supertype for all mass-metallicity relations. """
abstract type AbstractMZR{T <: Real} end
Base.Broadcast.broadcastable(m::AbstractMZR) = Ref(m)
(mzr::AbstractMZR)(Mstar::Real) = OH(mzr, Mstar) # Make all subtypes callable

"""
    OH(mzr::AbstractMZR, Mstar::Real)

Returns ``12 + \\log(\\text{O} / \\text{H})`` from the mass-metallicity relation `mzr` at stellar mass `Mstar` [M⊙].

# Examples
```jldoctest
julia> OH(AstroScalingRelations.Andrews2013, 1e8) ≈ 8.119245564173369
true
```
"""
OH(mzr::AbstractMZR, Mstar::Real)
OH(mzr::AbstractMZR, Mstar::u.Mass) = OH(mzr, u.ustrip(ua.Msun, Mstar))
"""
    dOH_dMstar(mzr::AbstractMZR, Mstar::Union{Real, Unitful.Mass})

Returns the partial derivative of the mass-metallicity relation `mzr` with respect to the stellar mass evaluated at `Mstar` [M⊙], where the metallicity is defined as ``12 + \\log(\\text{O} / \\text{H})``.

# Examples
```jldoctest
julia> dOH_dMstar(AstroScalingRelations.Andrews2013, 1e8) ≈ 2.1971001283624613e-9
true
```
"""
dOH_dMstar(mzr::AbstractMZR, Mstar::Real)
dOH_dMstar(mzr::AbstractMZR, Mstar::u.Mass) = dOH_dMstar(mzr, u.ustrip(ua.Msun, Mstar))

"""
    dOH_dlnMstar(mzr::AbstractMZR, Mstar::Union{Real, Unitful.Mass})

Returns the partial derivative of the mass-metallicity relation `mzr` with respect to the natural logarithm of the stellar mass evaluated at `Mstar` [M⊙], where the metallicity is defined as ``12 + \\log(\\text{O} / \\text{H})``. Mathematically, the following are equivalent: `dOH_dlnMstar(mzr, Mstar) = Mstar * dOH_dMstar(mzr, Mstar)`.

# Examples
```jldoctest
julia> dOH_dlnMstar(AstroScalingRelations.Andrews2013, 1e8) ≈ 2.1971001283624613e-9 * 1e8
true
```
"""
dOH_dlnMstar(mzr::AbstractMZR, Mstar::Real) = dOH_dMstar(mzr, Mstar) * Mstar
dOH_dlnMstar(mzr::AbstractMZR, Mstar::u.Mass) = dOH_dMstar(mzr, Mstar) * u.ustrip(ua.Msun, Mstar)

"""
    dOH_dlog10Mstar(mzr::AbstractMZR, Mstar::Union{Real, Unitful.Mass})

Returns the partial derivative of the mass-metallicity relation `mzr` with respect to the base-10 logarithm of the stellar mass evaluated at `Mstar` [M⊙], where the metallicity is defined as ``12 + \\log(\\text{O} / \\text{H})``. Mathematically, the following are equivalent: `dOH_dlog10Mstar(mzr, Mstar) = log(10) * dOH_dlnMstar(mzr, Mstar)`.

# Examples
```jldoctest
julia> dOH_dlog10Mstar(AstroScalingRelations.Andrews2013, 1e8) ≈ 2.1971001283624613e-9 *
                                                                 1e8 * log(10)
true
```
"""
dOH_dlog10Mstar(mzr::AbstractMZR, Mstar::Union{Real, u.Mass}) = dOH_dlnMstar(mzr, Mstar) * logten

"""
    OH_from_Z(Z::Real, X::Real, f_O::Real)

Returns ``12 + \\log(\\text{O} / \\text{H})`` calculated from the metal mass fraction `Z` assuming a hydrogen mass fraction `X` and that the fraction of metal mass that is oxygen is `f_O`.

```math
12 + \\log(\\text{O} / \\text{H}) = 12 + \\log \\left( \\frac{f_\\text{O} \\, Z}{16 X} \\right)
```

# Examples
```jldoctest
julia> OH_from_Z(0.0127, 0.74, 0.35) ≈ 8.574520062919332
true
```
"""
OH_from_Z(Z::Real, X::Real, f_O::Real) = 12 + log10(f_O * Z / 16X)

"""
    Z_from_OH(OH::Real, X::Real, f_O::Real)

Returns the metal mass fraction `Z` calculated from the oxygen abundance `OH` defined as ``12 + \\log(\\text{O} / \\text{H})`` assuming a hydrogen mass fraction `X` and that the fraction of metal mass that is oxygen is `f_O`.

```math
\\begin{aligned}
\\text{OH} &\\equiv 12 + \\log(\\text{O} / \\text{H}) \\newline
Z &= 10^{\\left( \\text{OH} - 12 \\right)} * 16 X / f_\\text{O}
\\end{aligned}
```

# Examples
```jldoctest
julia> 0.0127 ≈ Z_from_OH(OH_from_Z(0.0127, 0.74, 0.35),
                          0.74, 0.35)
true
```
"""
Z_from_OH(OH::Real, X::Real, f_O::Real) = exp10(OH - 12) * 16X / f_O

"""
    dZ_dOH(OH::Real, X::Real, f_O::Real)

Returns the partial derivative of the metal mass fraction `Z` with respect to the oxygen abundance `OH` defined as ``12 + \\log(\\text{O} / \\text{H})`` assuming a hydrogen mass fraction `X` and that the fraction of metal mass that is oxygen is `f_O`.

```math
\\begin{aligned}
\\text{OH} &\\equiv 12 + \\log(\\text{O} / \\text{H}) \\newline
Z &= 10^{\\left( \\text{OH} - 12 \\right)} * 16 X / f_\\text{O} \\newline
\\frac{\\partial \\, Z}{\\partial \\, \\text{OH}} &= \\ln(10) \\, 10^{\\left( \\text{OH} - 12 \\right)} * 16 X / f_\\text{O}
\\end{aligned}
```

# Examples
```jldoctest
julia> dZ_dOH(8.574520062919332, 0.74, 0.35) ≈ 0.029242830681024394
true
```
"""
dZ_dOH(OH::Real, X::Real, f_O::Real) = logten * exp10(OH - 12) * 16X / f_O

"""
    dZ_dMstar(mzr::AbstractMZR, Mstar::Union{Real, u.Mass}, X::Real, f_O::Real)

Returns the partial derivative of the metal mass fraction `Z` with respect to galaxy stellar mass `Mstar` [M⊙] assuming the mass-metallicity relation `mzr`, a hydrogen mass fraction `X`, and that the fraction of metal mass that is oxygen is `f_O`. Calculated as `oh = OH(mzr, Mstar); dZ_dMstar = dZ_dOH(oh, X, f_O) * dOH_dMstar(mzr, Mstar)`.

# Examples
```jldoctest
julia> dZ_dMstar(AstroScalingRelations.Andrews2013, 1e7, 0.74, 0.35) ≈ 7.33869576597908e-11
true
```
"""
function dZ_dMstar(mzr::AbstractMZR, Mstar::Union{Real, u.Mass}, X::Real, f_O::Real)
    oh = OH(mzr, Mstar)
    return dZ_dOH(oh, X, f_O) * dOH_dMstar(mzr, Mstar)
end
# using ForwardDiff: derivative; using AstroScalingRelations: Z_from_OH, OH, Andrews2013; derivative(mstar -> Z_from_OH(OH(Andrews2013, mstar), 0.74, 0.35), 1e7)


#### Types

"""
    MoustakasMZR(OH0::Real, M_to::Real, γ::Real)

General form for the asymptotic logarithmic MZR suggested by [Moustakas et al. 2011](https://arxiv.org/abs/1112.3300), given by

```math
12 + \\log(\\text{O} / \\text{H}) = \\text{OH0} - \\log \\left( 1 + \\left( \\frac{M_\\text{to}}{M_*} \\right)^\\gamma \\right)
```

where `OH0` is the asymptotic metallicity, `M_to` is the turnover mass, and `γ` is the power-law slope of the MZR at low stellar masses.

# Specific Instances
This MZR has been used in several publications that fit the parameters to data. These published results are available as the following:
 - `Andrews2013`: Result from [Andrews & Martini 2013](https://ui.adsabs.harvard.edu/abs/2013ApJ...765..140A/abstract) with `OH0 = 8.798`, `log10(M_to) = 8.901`, and `γ=0.640` based on fitting data from SDSS DR7 for galaxies with `log10(Mstar)` between 7.4 -- 10.5.
"""
struct MoustakasMZR{T <: Real} <: AbstractMZR{T}
    OH0::T
    M_to::T
    γ::T
end
MoustakasMZR(OH0::Real, M_to::Real, γ::Real) =
    MoustakasMZR(promote(OH0, M_to, γ)...)

OH(mzr::MoustakasMZR, Mstar::Real) = mzr.OH0 - log10(1 + (mzr.M_to/Mstar)^mzr.γ)
function dOH_dMstar(mzr::MoustakasMZR, Mstar::Real)
    OH0, M_to, γ = mzr.OH0, mzr.M_to, mzr.γ
    pl = (M_to/Mstar)^γ
    return  pl * γ / ((Mstar * logten) * (1 + pl))
end

"Result from [Andrews & Martini 2013](https://ui.adsabs.harvard.edu/abs/2013ApJ...765..140A/abstract) with `OH0 = 8.798`, `log10(M_to) = 8.901`, and `γ=0.640` based on fitting data from SDSS DR7 for galaxies with `log10(Mstar)` between 7.4 -- 10.5."
const Andrews2013 = MoustakasMZR(8.798, exp10(8.901), 0.640)

#################################################################################

"""
    CurtiMZR(OH0::Real, M_to::Real, γ::Real, β::Real)

General form for the asymptotic logarithmic MZR with adjustable saturation transition suggested by [Curti et al. 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.491..944C/abstract), given by

```math
12 + \\log(\\text{O} / \\text{H}) = \\text{OH0} - \\frac{\\gamma}{\\beta} \\, \\log \\left( 1 + \\left( \\frac{M_*}{M_\\text{to}} \\right)^{-\\beta} \\right)
```

Relative to the [`MoustakasMZR`](@ref), this MZR has one more free parameter which adjusts how broad the transition region is from the power-law MZR at low stellar masses to the asymptotic MZR at high stellar masses.

# Specific Instances
This MZR has been used in several publications that fit the parameters to data. These published results are available as the following:
 - `Curti2020`: Result from the paper that introduces this functional form with `OH0=8.793 ± 0.005`, `log10(M_to) = 10.02 ± 0.09`, `γ = 0.28 ± 0.02`, and `β = 1.2 ± 0.2`.
"""
struct CurtiMZR{T <: Real} <: AbstractMZR{T}
    OH0::T
    M_to::T
    γ::T
    β::T
end
CurtiMZR(OH0::Real, M_to::Real, γ::Real, β::Real) =
    CurtiMZR(promote(OH0, M_to, γ, β)...)

OH(mzr::CurtiMZR, Mstar::Real) = mzr.OH0 - mzr.γ / mzr.β * log10(1 + (Mstar/mzr.M_to)^-mzr.β)
function dOH_dMstar(mzr::CurtiMZR, Mstar::Real)
    OH0, M_to, γ, β = mzr.OH0, mzr.M_to, mzr.γ, mzr.β
    r1 = Mstar/M_to
    r2 = r1^-β
    return r2/r1 * γ / ((M_to * logten) * (1 + r2))
end

"Result from [Curti et al. 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.491..944C/abstract) with `OH0=8.793 ± 0.005`, `log10(M_to) = 10.02 ± 0.09`, `γ = 0.28 ± 0.02`, and `β = 1.2 ± 0.2` based on fitting data from SDSS DR7 -- this is mostly the same data as was used to derive the `Andrews2013` result, but Curti 2020 used new strong-line metallicity calibrations."
const Curti2020 = CurtiMZR(8.793, exp10(10.02), 0.28, 1.2)

#################################################################################

# struct Ma2016MZR{T <: Real} <: AbstractMZR{T}
#     γ::T    # Power law slope
#     # Zten::T # Metal mass fraction at M_* = 10^10 solar masses
#     A::T    # A, B, C control redshift evolution of Z_{10}, which is
#     B::T    # the metal mass fraction at M_* = 10^10 solar masses
#     C::T
# end
# Ma2016MZR(γ::Real, A::Real, B::Real, C::Real) =
#     Ma2016MZR(promote(γ, A, B, C)...)
"""
    PowerLawMZR(α::Real, β::Real, logMstar0::Real)
Mass-metallicity model described by a single power law index `α`, a metallicity normalization `β`, and the logarithm of the stellar mass `logMstar0 = log10(Mstar0 [M⊙])` at which the metallicity is `β`. Such a power law MZR is often used when extrapolating literature results to low masses, e.g., ``\\text{M}_* < 10^8 \\; \\text{M}_\\odot.``

# Specific Instances
 - [`Ma2016Gas`](@ref) and [`Ma2016Stars`](@ref) are redshift-dependent power law MZR models fit to the gas-phase and stellar mass-metallicity relations in simulated FIRE galaxies.
"""
struct PowerLawMZR{T <: Real} <: AbstractMZR{T}
    α::T # Power law slope
    β::T # Metallicity at log10(Mstar) = logMstar0
    logMstar0::T 
end
PowerLawMZR(α::Real, β::Real, logMstar0::Real) = 
    PowerLawMZR(promote(α, β, logMstar0)...)
OH(mzr::PowerLawMZR, Mstar::Real) = mzr.α * (log10(Mstar) - mzr.logMstar0) + mzr.β
dOH_dMstar(mzr::PowerLawMZR, Mstar::Real) = mzr.α / Mstar / logten

"""
    Ma2016Gas(z::Real, γ::Real=0.35; A::Real=0.93, B::Real=0.43, C::Real=-1.05)::PowerLawMZR
Returns a [`PowerLawMZR`](@ref) instance valid at redshift `z` given the gas-phase mass-metallicity relation measured for FIRE galaxies by [Ma et al. 2016](https://ui.adsabs.harvard.edu/abs/2016MNRAS.456.2140M/abstract). `γ` is the power law slope and `A, B, C` control the redshift evolution of the model.

See also [`Ma2016Stars`](@ref) for the stellar mass-metallicity relation measured in the same work. 
"""
function Ma2016Gas(z::Real, γ::Real=0.35; A::Real=0.93, B::Real=0.43, C::Real=-1.05)
    Zten = A * exp(-B * z) + C
    # Zten is the log(Zgas / Z⊙) = 12 + log(O/H) - 9; so we have to add 9
    return PowerLawMZR(γ, Zten + 9, 10)
end
"""
    Ma2016Stars(z::Real, γ::Real=0.4; A::Real=0.67, B::Real=0.5, C::Real=-1.04)::PowerLawMZR
Returns a [`PowerLawMZR`](@ref) instance valid at redshift `z` given the stellar mass-metallicity relation measured for FIRE galaxies by [Ma et al. 2016](https://ui.adsabs.harvard.edu/abs/2016MNRAS.456.2140M/abstract). `γ` is the power law slope and `A, B, C` control the redshift evolution of the model.

See also [`Ma2016Gas`](@ref) for the gas-phase mass-metallicity relation measured in the same work. 
"""
function Ma2016Stars(z::Real, γ::Real=0.4; A::Real=0.67, B::Real=0.5, C::Real=-1.04)
    Zten = A * exp(-B * z) + C
    # Zten is the log(Zgas / Z⊙) = 12 + log(O/H) - 9; so we have to add 9
    return PowerLawMZR(γ, Zten + 9, 10)
end
