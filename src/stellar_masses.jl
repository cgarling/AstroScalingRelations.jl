# Functions for estimating stellar masses of galaxies from photometric observables.

"""
    Longhetti2009(mK, VK, z, d_mod)

Returns the stellar mass in solar masses for an early-type galaxy using the empirical
near-infrared mass estimator from Section 5 of
[Longhetti & Saracco 2009](https://doi.org/10.1111/j.1365-2966.2008.14375.x).

The estimator requires only the apparent K-band magnitude and the observed (V-K) colour,
making it independent of multi-wavelength SED fitting. It was calibrated on early-type
(elliptical) galaxies in the redshift range ``1 < z < 2``.

# Arguments
- `mK`: apparent K-band magnitude of the galaxy
- `VK`: observed (V-K) colour of the galaxy
- `z`: redshift of the galaxy (valid range: ``1 < z < 2``)
- `d_mod`: distance modulus in magnitudes (cosmology-dependent; e.g., ``d_\\text{mod} = 5 \\log_{10}(d_L / 10\\,\\text{pc})``)

# Method
The estimator combines Equations 14 and 15 from the paper:

```math
\\begin{aligned}
\\log_{10}(M/L_K) &= 0.086 \\times (V-K)_\\text{obs} - 0.863 \\\\
M_K &= m_K - d_\\text{mod}(z) - k_\\text{cor}(z) \\\\
k_\\text{cor}(z) &= -0.32 - 0.65\\,z - 0.314\\,z^2 - 0.0735\\,z^3 \\\\
\\log_{10}(M/M_\\odot) &= \\log_{10}(M/L_K) - 0.4 \\, M_K + 1.364
\\end{aligned}
```

where ``0.4 \\times M^\\odot_K = 1.364`` (with ``M^\\odot_K = 3.41``) and the k-correction
polynomial uses the coefficients from Table 3, case i) (all models except CB08 and Ma05,
at solar metallicity).

# Notes
- The scatter of the ``\\log_{10}(M/L_K)`` vs. (V-K) relation (Eq. 15) is ``\\sigma = 0.04``
  with a maximum deviation of ``\\pm 0.1``.
- Overall stellar mass estimates are recovered to within a factor of ``\\sim 0.8``--``1.1``
  for 68% of simulated early-type galaxies at ``1 < z < 2``.
- The Salpeter IMF was used when calibrating the M/L models underlying this estimator.

# Examples
```jldoctest
julia> Longhetti2009(20.0, 5.0, 1.5, 45.0) ≈ 1.0744222661470343e10
true

julia> Longhetti2009(19.0, 4.5, 1.0, 44.1) ≈ 2.4266100950824142e10
true
```
"""
function Longhetti2009(mK, VK, z, d_mod)
    kcor = -0.32 - 0.65*z - 0.314*z^2 - 0.0735*z^3
    MK = mK - d_mod - kcor
    log_ML = 0.086 * VK - 0.863
    log_M = log_ML - 0.4 * MK + 1.364
    return exp10(log_M)
end
