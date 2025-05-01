# ThermovisorData

[![Build Status](https://github.com/Manarom/ThermovisorData.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Manarom/ThermovisorData.jl/actions/workflows/CI.yml?query=branch%3Amaster)

**ThervisorData.jl** is a small package designed to process static thermal images stored as matrices in CSV format, where each matrix element represents a temperature value. It enables users to load images from files, calculate temperature distributions, and compute statistical analyses for temperatures along specified lines. It also calculates averaged angular and radial temperature distributions (along with standard deviations) within Regions of Interest (ROIs) such as circles, squares, and rectangles. These ROI objects can be fitted to thermally distinct areas (relative to their surroundings), such as the most heated regions within the scene.

 In `/test` folder there are two [Pluto](https://plutojl.org/) notebooks, which can be used as examples of package usage.


  Full documentation is available at https://manarom.github.io/BandPyrometry.jl/

# Installation 

1) Download [Julia](https://julialang.org/downloads)

#### For usage

2a) Clone this repository to your local machine 

3a) use `include("BandPyrometry_folder\src\BandPyrometry.jl)` in REPL or other module's to bring BandPyrometry module to the global scope

4a) use `using .BandPyrometry` to bring the module and its content to the corresponding namespace

#### For development

2c) Clone this repository to `username/.julia/dev/`.

3c) Enter the package manager in REPL by pressing `]`  then add the package by typing `dev BandPyrometry`
