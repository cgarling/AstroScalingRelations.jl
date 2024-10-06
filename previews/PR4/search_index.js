var documenterSearchIndex = {"docs":
[{"location":"public_methods/#methods","page":"Public Methods","title":"Public Methods","text":"","category":"section"},{"location":"public_methods/","page":"Public Methods","title":"Public Methods","text":"The following methods are part of our publicly exported API.","category":"page"},{"location":"public_methods/#Stellar-Mass-Halo-Mass-Relation","page":"Public Methods","title":"Stellar-Mass-Halo-Mass Relation","text":"","category":"section"},{"location":"public_methods/","page":"Public Methods","title":"Public Methods","text":"flum\nMoster2013\nBehroozi2013\nGK14\nGK17_field\nGK17_sat","category":"page"},{"location":"public_methods/#AstroScalingRelations.flum","page":"Public Methods","title":"AstroScalingRelations.flum","text":"flum(Mh::Real, mdef::String=\"vir\")\n\nAnalytic approximations to the fraction of halos with mass definition mdef of mass Mh (in solar masses) that are luminous, based on Fig. 3 and 12 of Dooley 2017b. If mdef=\"late z\", then it is most similar to the purple dashed or solid green models in Figure 12. This was Annika's code. mdef=\"200c\" is most similar to the Barber (blue) model in Figure 3. mdef=\"vir\" is most similar to the yellow model in Figure 3. mdef=\"350c\" is most similar to the red model in Figure 3. If mdef=\"const\", then just return 1; e.g. 100% luminous fraction at all halo masses.\n\nExamples\n\njulia> flum(1e10, \"late z\") ≈ 0.9999938558253978\ntrue\n\njulia> flum(1e10, \"200c\") ≈ 0.9994793959258891\ntrue\n\njulia> flum(1e10, \"vir\") ≈ 0.9878715650157257\ntrue\n\njulia> flum(1e10, \"350c\") ≈ 0.9933071490757152\ntrue\n\njulia> flum(1e10, \"const\") ≈ 1.0\ntrue\n\n\n\n\n\n","category":"function"},{"location":"public_methods/#AstroScalingRelations.Moster2013","page":"Public Methods","title":"AstroScalingRelations.Moster2013","text":"Moster2013(Mh, z, M10=11.59, M11=1.195, N10=0.0351, N11=-0.0247, B10=1.376, B11=-0.826, gamma10=0.608, gamma11=0.329)\n\nGiven a halo mass Mh (M200c in solar masses) and a redshift z, return the median stellar mass from the Moster et al. 2013 empirical model. If Mh is a Unitful.Mass, will return a Unitful.Mass. \n\nExamples\n\njulia> Moster2013( 1e12, 0.0 ) ≈ 3.427515185606664e10\ntrue\n\njulia> Moster2013( 1e12 * UnitfulAstro.Msun, 0.0 ) ≈ 3.427515185606664e10 * UnitfulAstro.Msun\ntrue\n\nNotes\n\nMoster et al. 2013 utilized the Millennium simulation (Springel et al. 2005) for high-mass halos and the Millennium-II simulation (Boylan-Kolchin et al. 2009) for low-mass halos. At low redshift, they use the stellar mass function of Li & White 2009 derived from SDSS DR7, ranging from stellar masses of 10^8.3 to 10^11.7 solar masses. This is supplemented by the low-mass stellar mass function from Baldry et al. 2008 for stellar masses between 10^7.4 to 10^8.3 solar masses. The high-redshift galaxy sample covers redshifts 0 < z < 4, using data from Perez-Gonzalez et al. 2008 for high-mass galaxies and Santini et al. 2012 for low-mass galaxies. See Figure 6 in Moster et al. 2013 for a visualization of the adopted stellar mass functions. \n\n\n\n\n\n","category":"function"},{"location":"public_methods/#AstroScalingRelations.Behroozi2013","page":"Public Methods","title":"AstroScalingRelations.Behroozi2013","text":"Behroozi2013(Mh, z)\n\nGiven a halo mass Mh (specifically the peak virial halo mass in solar masses, as defined in Bryan & Norman 1998) and a redshift z, return the median stellar mass from the Behroozi et al. 2013 empirical model. If Mh is a Unitful.Mass, will return a Unitful.Mass.\n\nExamples\n\njulia> Behroozi2013( 1e12, 0.0 ) ≈ 2.6798246456860065e10\ntrue\n\njulia> Behroozi2013( 1e12 * UnitfulAstro.Msun, 0.0 ) ≈ 2.6798246456860065e10 * UnitfulAstro.Msun\ntrue\n\n\n\n\n\n","category":"function"},{"location":"public_methods/#AstroScalingRelations.GK14","page":"Public Methods","title":"AstroScalingRelations.GK14","text":"GK14(Mh, z)\n\nGiven a halo mass Mh (specifically the peak virial halo mass in solar masses, as defined in Bryan & Norman 1998) and a redshift z, return the median stellar mass from the Garrison-Kimmel et al. 2014 empirical model. If Mh is a Unitful.Mass, will return a Unitful.Mass. This is the same model as Behroozi2013, but with a shallower low-mass slope (α). See Figures 9 and 10. \n\nExamples\n\njulia> GK14( 1e12, 0.0 ) ≈ 2.894144183405201e10\ntrue\n\njulia> GK14( 1e12 * UnitfulAstro.Msun, 0.0 ) ≈ 2.894144183405201e10 * UnitfulAstro.Msun\ntrue\n\n\n\n\n\n","category":"function"},{"location":"public_methods/#AstroScalingRelations.GK17_field","page":"Public Methods","title":"AstroScalingRelations.GK17_field","text":"GK17_field(Mh, z, nu=-0.2)\n\nGiven a halo mass Mh (specifically the peak virial halo mass in solar masses, as defined in Bryan & Norman 1998) and a redshift z, return the median stellar mass for field galaxies from the Garrison-Kimmel et al. 2017 empirical model. If Mh is a Unitful.Mass, will return a Unitful.Mass. This uses the model of Behroozi2013 at high masses, but has a variable low-mass slope (α) and scatter (σν). The argument nu is the slope of the scatter in the SMHM relation (ν in their Equation 3). The default is -0.2 as in Dooley et al. 2017 \"An Observers Guide ...\".\n\nExamples\n\njulia> GK17_field( 1e12, 0.0, -0.2 ) ≈ 2.9300795754861504e10\ntrue\n\njulia> GK17_field( 1e12 * UnitfulAstro.Msun, 0.0, -0.2 ) ≈ 2.9300795754861504e10 * UnitfulAstro.Msun\ntrue\n\n\n\n\n\n","category":"function"},{"location":"public_methods/#AstroScalingRelations.GK17_sat","page":"Public Methods","title":"AstroScalingRelations.GK17_sat","text":"GK17_sat(Mh, z, nu=-0.2)\n\nGiven a halo mass Mh (specifically the peak virial halo mass in solar masses, as defined in Bryan & Norman 1998) and a redshift z, return the median stellar mass for satellite galaxies from the Garrison-Kimmel et al. 2017 empirical model. If Mh is a Unitful.Mass, will return a Unitful.Mass. This uses the model of Behroozi2013 at high masses, but has a variable low-mass slope (α) and scatter (σν). The argument nu is the slope of the scatter in the SMHM relation (ν in their Equation 3). The default is -0.2 as in Dooley et al. 2017 \"An Observers Guide ...\".\n\nExamples\n\njulia> GK17_sat( 1e12, 0.0, -0.2 ) ≈ 2.860110341918958e10\ntrue\n\njulia> GK17_sat( 1e12 * UnitfulAstro.Msun, 0.0, -0.2 ) ≈ 2.860110341918958e10 * UnitfulAstro.Msun\ntrue\n\n\n\n\n\n","category":"function"},{"location":"public_methods/#Gas-Outflow-Models-(Mass-Loading-Factors)","page":"Public Methods","title":"Gas Outflow Models (Mass-Loading Factors)","text":"","category":"section"},{"location":"public_methods/","page":"Public Methods","title":"Public Methods","text":"Muratov2015\nChristensen2016\nPandya2021","category":"page"},{"location":"public_methods/#AstroScalingRelations.Muratov2015","page":"Public Methods","title":"AstroScalingRelations.Muratov2015","text":"Muratov2015(Mstar)\nMuratov2015(Vvir, z)\n\nThe model for the galaxy-scale mass-loading factor η from Muratov et al. 2015, based on FIRE-1 simulations. Two parameterizations are allowed; one which is based on the galaxy's stellar mass (Mstar) in solar masses and one which is based on the circular velocity at the virial radius (Vvir) of the galaxy's dark matter halo in km/s and the redshift of evaluation (z). Both parameterizations can be Unitful quantities. See also Pandya2021. \n\nExamples\n\njulia> Muratov2015( 1e10 ) ≈ 3.6\ntrue\n\njulia> Muratov2015( 1e10 * UnitfulAstro.Msun ) ≈ 3.6\ntrue\n\njulia> Muratov2015( 40, 1 ) ≈ 26.135392171726632\ntrue\n\njulia> Muratov2015( 40 * Unitful.km / Unitful.s, 1 ) ≈ 26.135392171726632\ntrue\n\n\n\n\n\n","category":"function"},{"location":"public_methods/#AstroScalingRelations.Christensen2016","page":"Public Methods","title":"AstroScalingRelations.Christensen2016","text":"Christensen2016(Vvir)\n\nThe model for the halo-scale mass-loading factor η from Christensen et al. 2016, based on simulations run with the GASOLINE hydrodynamical model. The model requires the circular velocity at the virial radius (Vvir) of the galaxy's dark matter halo in km/s. Vvir can be a Unitful.Velocity. Note that the normalization of the relation is not given explicitly in the paper, only the scaling is; the prefactor here was estimated from their Figure 11. \n\nExamples\n\njulia> Christensen2016( 40 ) ≈ 3.0028111668298063\ntrue\n\njulia> Christensen2016( 40 * Unitful.km / Unitful.s ) ≈ 3.0028111668298063\ntrue\n\n\n\n\n\n","category":"function"},{"location":"public_methods/#AstroScalingRelations.Pandya2021","page":"Public Methods","title":"AstroScalingRelations.Pandya2021","text":"Pandya2021(Mstar)\nPandya2021(Vvir, z)\n\nThe model for the galaxy-scale mass-loading factor η from Pandya et al. 2021, based on FIRE-2 simulations. Two parameterizations are allowed; one which is based on the galaxy's stellar mass (Mstar) in solar masses and one which is based on the circular velocity at the virial radius (Vvir) of the galaxy's dark matter halo in km/s and the redshift of evaluation (z). Originally calculated for redshift 0  z  4 but will extrapolate outside this range. Both parameterizations can be Unitful quantities. This work finds that halo-scale mass-loading factors are larger than galaxy-scale mass-loading factors due to entrainment of additional CGM gas (see Section 5 and Figure 14). \n\nExamples\n\njulia> isapprox( Pandya2021( 1e10 ), 0.6309573444801929; rtol=1e-7 )\ntrue\n\njulia> isapprox( Pandya2021( 1e10 * UnitfulAstro.Msun ), 0.6309573444801929; rtol=1e-7 )\ntrue\n\njulia> all( isapprox.( Pandya2021.( 40, (0.0, 1.0, 3.0) ), (12.977828436995544, 21.609465380663718, 24.88169815959356); rtol=1e-7 ) )\ntrue\n\njulia> isapprox( Pandya2021( 40 * Unitful.km / Unitful.s, 1 ), 21.609465380663718; rtol=1e-7 )\ntrue\n\n\n\n\n\n","category":"function"},{"location":"public_methods/#Gas-Mass-/-Profile-Relations","page":"Public Methods","title":"Gas Mass / Profile Relations","text":"","category":"section"},{"location":"public_methods/","page":"Public Methods","title":"Public Methods","text":"Papastergis2012\nBradford2015\nScoville2017\nWang2016_DHI\nWang2016_rho_rs","category":"page"},{"location":"public_methods/#AstroScalingRelations.Papastergis2012","page":"Public Methods","title":"AstroScalingRelations.Papastergis2012","text":"Papastergis2012(Mstar)\n\nReturns the mean neutral hydrogen gas mass of a galaxy with stellar mass Mstar in solar masses; this is the black line in Figure 19 of Papastergis et al. 2012. Based on a crossmatch between SDSS galaxies and sources detected in the blind ALFALFA α.40 survey. If Mstar is a Unitful.Mass, will return a Unitful.Mass.\n\nExamples\n\njulia> Papastergis2012(1e8) ≈ 2.0417379446695232e8\ntrue\n\njulia> Papastergis2012(1e8 * UnitfulAstro.Msun) ≈ 2.0417379446695232e8 * UnitfulAstro.Msun\ntrue\n\n\n\n\n\n","category":"function"},{"location":"public_methods/#AstroScalingRelations.Bradford2015","page":"Public Methods","title":"AstroScalingRelations.Bradford2015","text":"Bradford2015(Mstar)\n\nReturns the median neutral hydrogen gas mass of an isolated galaxy with stellar mass Mstar in solar masses according to Equation in Bradford et al. 2015. Uses coefficients from the first two rows of Table 3. These are derived for galaxies selected from SDSS DR8 with stellar masses between exp10(7) and exp10(9.5).\n\nExamples\n\njulia> Bradford2015(1e8) ≈ 4.487453899331332e8\ntrue\n\njulia> Bradford2015(1e8 * UnitfulAstro.Msun) ≈ 4.487453899331332e8 * UnitfulAstro.Msun\ntrue\n\n\n\n\n\n","category":"function"},{"location":"public_methods/#AstroScalingRelations.Scoville2017","page":"Public Methods","title":"AstroScalingRelations.Scoville2017","text":"Scoville2017(Mstar [, z])\n\nReturns the median neutral hydrogen gas mass of a galaxy with stellar mass Mstar in solar masses at redshift z according to Equation 6 in Scoville et al. 2017. This does not include the specfic star formation rate dependence. These are derived from ALMA dust continuum emission (see Section 4) for well-studied galaxies in the COSMOS field between redshifts 0.3 and 4.5. If Mstar is a Unitful.Mass, will return a Unitful.Mass.\n\nExamples\n\njulia> Scoville2017(1e8) ≈ 1.7759037070772731e9\ntrue\n\njulia> Scoville2017(1e8 * UnitfulAstro.Msun) ≈  1.7759037070772731e9 * UnitfulAstro.Msun\ntrue\n\njulia> Scoville2017(1e8, 1.0) ≈ 6.357913365552342e9\ntrue\n\njulia> Scoville2017(1e8 * UnitfulAstro.Msun, 1.0) ≈ 6.357913365552342e9 * UnitfulAstro.Msun\ntrue\n\n\n\n\n\n","category":"function"},{"location":"public_methods/#AstroScalingRelations.Wang2016_DHI","page":"Public Methods","title":"AstroScalingRelations.Wang2016_DHI","text":"Wang2016_DHI(MHI)\n\nReturns the HI diameter, defined as the diameter at which the HI surface density equals 1 Msun/pc^2, in kiloparsecs given an HI mass follwing Equation 2 in Wang et al. 2016. If Mstar is a Unitful.Mass, will return a Unitful.Length.\n\nExamples\n\njulia> Wang2016_DHI(1e8) ≈ 5.688529308438413\ntrue\n\njulia> Wang2016_DHI(1e8 * UnitfulAstro.Msun) ≈ 5.688529308438413 * UnitfulAstro.kpc\ntrue\n\n\n\n\n\n","category":"function"},{"location":"public_methods/#AstroScalingRelations.Wang2016_rho_rs","page":"Public Methods","title":"AstroScalingRelations.Wang2016_rho_rs","text":"(Σ₀, rs) = Wang2016_rho_rs(MHI)\n\nReturns the central surface density Σ₀ in Msun/pc^2 and the scale radius rs in parsecs that defines the surface density profile Σ(r) = Σ₀ exp (-rrs) for a galaxy with an HI disk mass of MHI in solar masses, using the Wang et al. 2016 relation between HI mass and HI disk diameter.\n\nExamples\n\njulia> all( Wang2016_rho_rs(1e6) .≈ (4.009796759878976, 199.2273166291845) )\ntrue\n\njulia> all( Wang2016_rho_rs(1e12) .≈ (7.045959041825033, 150293.43437081948) )\ntrue\n\n\n\n\n\n","category":"function"},{"location":"public_methods/#Galaxy-Sizes","page":"Public Methods","title":"Galaxy Sizes","text":"","category":"section"},{"location":"public_methods/","page":"Public Methods","title":"Public Methods","text":"galaxy_size","category":"page"},{"location":"public_methods/#AstroScalingRelations.galaxy_size","page":"Public Methods","title":"AstroScalingRelations.galaxy_size","text":"galaxy_size(Rvir; Aᵣ=0.02, fₚfₖ=0.78) = Aᵣ * fₚfₖ * Rvir\ngalaxy_size(Mh, ρthresh; Aᵣ=0.02, fₚfₖ=0.78) = Aᵣ * fₚfₖ * cbrt( 3 * Mh / (4π * ρthresh) )\n\nScaling relation for galaxy half-light radius with the virial radius of the dark matter halo such that\n\nbeginaligned\nρ_textthresh = ρ_c(z)  Δ_textvir newline\nR_textvir = sqrt3 frac3  textM_h4π  ρ_c(z)   Δ_textvir newline\nR_h(textM_h z) = A_r  f_p f_k  R_textvir\nendaligned\n\nThe units of the returned R_h are the same as the provided R_textvir, or if using the (Mh, ρthresh) signature, the unit of cbrt( 3 * Mh / (4π * ρthresh) ).\n\nDefault keyword arguments Aᵣ=002 from Jiang et al. 2019 and fₚfₖ=078 for the projection correction for spheroidal systems from Somerville et al. 2018.\n\nExamples\n\njulia> galaxy_size(1e12, 12000; Aᵣ=0.02, fₚfₖ=0.78) ≈ 4.227023347086455\ntrue\n\n\n\n\n\n","category":"function"},{"location":"#guide","page":"AstroScalingRelations.jl Overview","title":"AstroScalingRelations.jl Overview","text":"","category":"section"},{"location":"","page":"AstroScalingRelations.jl Overview","title":"AstroScalingRelations.jl Overview","text":"This package implements astrophysical scaling relations like the stellar-mass-halo-mass function and stellar-mass-HI-mass relation. Methods are generally named after the paper they come from. Links to the relevant papers are included in the docstrings. See further documentation on the methods here. ","category":"page"}]
}