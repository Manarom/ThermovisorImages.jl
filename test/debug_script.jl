source_dir = joinpath(abspath(joinpath(@__DIR__, "..")),"src")

using Revise
includet(joinpath(@__DIR__,"ThermovisorData.jl"))

#include(joinpath(source_dir,"ThermovisorData.jl"))

# search all files in dir 
files_in_dir = ThermovisorData.find_temperature_files()

# selecting tag #1
selected_tag = collect(keys(files_in_dir))[1]

# loading image data (returns tuple of rescaled image and file creation date)
(img,date) = ThermovisorData.read_temperature_file(files_in_dir[selected_tag][2])
# creating `CentredObj` of various shapes
circle = ThermovisorData.CircleObj() # empty obj
# creating CentredObj using universal constructor
square = ThermovisorData.obj_from_vect(ThermovisorData.SquareObj,[200,100,100])
rect = ThermovisorData.obj_from_vect(ThermovisorData.RectangleObj,[50,80,100,100])
# simple wrapper around the image labeling functions from 
markers= ThermovisorData.marker_image(img)
flt = ThermovisorData.filter_image(img,markers)
ThermovisorData.fit_centred_obj!(circle,flt,nothing,fit_reduced=true)
ThermovisorData.draw!(copy(flt.reduced),circle)
ThermovisorData.draw!(copy(flt.reduced),square)
ThermovisorData.draw!(copy(flt.reduced),rect,fill=false)

using Images, ImageShow
ThermovisorData.diagonal_points(rect)
simshow(ThermovisorData.fill_im!(copy(flt.reduced),rect))

ThermovisorData.draw(circle)
circle
ThermovisorData.diagonal_points(circle)
im_copy = deepcopy(img)
ThermovisorData.along_line_distribution_xiaolin_wu(im_copy.initial, 1, 1, 200, 200)
simshow(im_copy.initial)


(along_line,distrib) = ThermovisorData.radial_distribution(flt.reduced,circle,0.0:180)
(L,D,stdD,l_b,u_b,t_values) = ThermovisorData.radial_distribution_statistics(along_line,distrib)
ThermovisorData.plot_radial_distribution_statistics(L,D,stdD,l_b,u_b)

(angs_p,Dang,stdDang,l_band,u_band,t_valuesang) = ThermovisorData.angular_distribution_statistics(collect(0.0:180),
            along_line,distrib)

ThermovisorData.plot_angular_distribution_statistics(angs_p,Dang,stdDang,l_band,u_band)


x0 = [0.0,0.0,0.0]
immm_flag = flt.reduced .>300.0
using ImageShow
simshow(immm_flag)
using BenchmarkTools
ThermovisorData.fill_x0!(x0,immm_flag,square)
ThermovisorData.fill_x0!(x0,immm_flag,circle)

ThermovisorData.fill_x0!([0.0,x0...],immm_flag,rect)