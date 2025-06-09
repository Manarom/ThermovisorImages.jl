push!(LOAD_PATH,"../src/")
#include("../src/ThermovisorImages.jl")
using Documenter, ThermovisorImages
mathengine = Documenter.MathJax3()
makedocs(
        sitename = "ThermovisorImages.jl",
        repo="https://github.com/Manarom/ThermovisorImages.jl/blob/{commit}{path}#{line}",
        highlightsig = false,
        checkdocs = :none,
        format=Documenter.HTML(size_threshold = nothing),
        pages=[
                "Home" => "index.md"
                "Examples"=>["Pluto notebooks"=>"pluto_tests_git.md"]
                "Modules" => ["ThermovisorImages" =>"thermovisordata.md"] 
               ]#
	)
deploydocs(;
        repo="github.com/Manarom/ThermovisorImages.jl", 
        devbranch = "master",
        devurl="dev",
        target = "build",
        branch = "gh-pages",
        versions = ["stable" => "v^", "v#.#" ]
    )