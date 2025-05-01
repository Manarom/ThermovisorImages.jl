push!(LOAD_PATH,"../src/")
include("../src/ThermovisorData.jl")
using Documenter,.ThermovisorData
mathengine = Documenter.MathJax3()
makedocs(
        sitename = "ThermovisorData.jl",
        highlightsig = false,
        checkdocs = :none,
        format=Documenter.HTML(size_threshold = nothing),
        pages=[
                "Home" => "index.md"
                #"Examples"=>["ThermovisorData-test"=>"thermovisordata_test.md"]
                "Modules" => ["ThermovisorData" =>"thermovisordata.md"] 
               ]#
	)
#deploydocs(;
#         repo="github.com/Manarom/BandPyrometry"
#)