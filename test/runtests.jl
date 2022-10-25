using AstroScalingRelations
using Test
using Documenter

DocMeta.setdocmeta!(AstroScalingRelations, :DocTestSetup, :(using AstroScalingRelations; import Unitful; import UnitfulAstro); recursive=true)
doctest(AstroScalingRelations,manual=false)
