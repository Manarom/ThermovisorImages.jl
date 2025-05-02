# ThermovisorData
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://manarom.github.io/ThermovisorData.jl)
[![Build Status](https://github.com/Manarom/ThermovisorData.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Manarom/ThermovisorData.jl/actions/workflows/CI.yml?query=branch%3Amaster)



**ThervisorData.jl** is a small package designed to process static thermal images stored as matrices. Each matrix element represents a temperature value. Package enables users to load images from csv-files, calculate temperature distributions, and compute statistical analyses for temperatures along specified lines. It also calculates averaged angular and radial temperature distributions (along with standard deviations) within Regions of Interest (ROIs) such as circles, squares, and rectangles. These ROI objects can be fitted to thermally distinct areas (relative to their surroundings), such as the most heated regions within the scene.

<img src="https://manarom.github.io/ThermovisorData.jl/master/docs/src/assets/initial_image.png" align="middle"  />


 In `/notebook` folder there is a [Pluto](https://plutojl.org/) notebook, which can be used as an example of package usage.


  Full documentation is available at  [documentation](https://manarom.github.io/ThermovisorData.jl/)

# Installation 

1) Download [Julia](https://julialang.org/downloads)

#### For usage

2a) Clone this repository to your local machine 

3a) use `include("thermovisor_data_folder\src\ThermovisorData.jl)` in REPL or other module's to bring BandPyrometry module to the global scope

4a) use `using .ThermovisorData` to bring the module and its content to the corresponding namespace

#### For development

2c) Clone this repository to `username/.julia/dev/`.

3c) Enter the package manager in REPL by pressing `]`  then add the package by typing `dev ThermovisorData`
