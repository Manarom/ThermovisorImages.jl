source_dir = joinpath(abspath(joinpath(@__DIR__, "..")),"src")

using Revise,ImageShow,Images,BenchmarkTools
includet(joinpath(source_dir,"ThermovisorData.jl"))
#=
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
ThermovisorData.fit_centred_obj!(circle,flt,fit_reduced=true)
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
=#

img = ThermovisorData.read_temperature_file(raw".\thermal images\2025-05-30 13-55-16_T600-2.csv")[1]

markers = ThermovisorData.marker_image(img)
simshow(markers)
mIm = ThermovisorData.MarkeredImage(markers)
#=
function sort_markers_by_area!(markers::Matrix{Int};total_number::Int = -1,rev::Bool=true)
    max_label = maximum(markers)

    if total_number <=0 || total_number >max_label 
        total_number = max_label
    end
    sums = Vector{Float64}(undef,max_label)
    inds = Vector{Int}(undef,max_label)
	flag = similar(markers,Bool)
    map!(i->i==1 ,flag, markers)
	v = @view markers[flag]
	ViewsVect = Vector{typeof(v)}(undef,total_number)
    ViewsVect[1] = v
    for i in 2:max_label
        map!(m->m==i ,flag, markers)
        ViewsVect[i] = @view markers[flag] 
    end
    map!(sum,sums,ViewsVect)
    #=for i in 1:max_label
        sums[i] = sum(ViewsVect[i])
    end=#
    sortperm!(inds,sums,rev=rev)
    for i = 1:total_number
        fill!(ViewsVect[inds[i]],i)
    end
    if total_number<max_label
        for i in total_number+1:max_label
            fill!(ViewsVect[inds[i]],0)
        end
    end
    return markers
end
=#


using Images, ImageSegmentation
use_test_image=true
if use_test_image
    img = load(download("http://docs.opencv.org/3.1.0/water_coins.jpg"));
    im_float =1.0 .-  @.  img |> Gray |>Float64
else

end
#im_float = 1 .-im_float
simshow(im_float)
maximum(im_float)
minimum(im_float)


thresh = otsu_threshold(1 .-im_float)
bw = im_float .> 0.5

dist = 1 .-distance_transform(feature_transform(bw));
maximum(dist)
minimum(dist)
simshow(dist)
markers = label_components(dist .< -8)
simshow(markers)
segments = watershed(dist, markers)

simshow( labels_map(segments) .* (0.0 .+ bw))       #shows segmented ima
simshow(markers)


im_rescaled = ThermovisorData.RescaledImage(im_float)
markers2 = ThermovisorData.marker_image(im_rescaled)
simshow(markers2)

simshow(segments.label_map)


##########################################
using Images, ImageSegmentation
img = load(download("http://docs.opencv.org/3.1.0/water_coins.jpg"));
im_float =1.0 .-  @.  img |> Gray |>Float64
simshow(im_float)
maximum(im_float)
minimum(im_float)

rescaled = ThermovisorData.RescaledImage(im_float)
marker = marker_image(rescaled)
simshow(marker)

function marker_image(rescaled::ThermovisorData.RescaledImage;
    level_threshold::Float64=-1.0,distance_threshold::Float64=-10.0)
if level_threshold>1 || level_threshold <=0 # if threshold is not settled explicitly the otsu thresholding algorithm is used
level_threshold = otsu_threshold(rescaled.im)
end
# thermal image is an image with several region of higher temperature
# we want to implement the watershed algorithm to separate patterns from each other
# first we need to negate image 
binarized =@. (1 - rescaled.im) > (1 -level_threshold)
dist = distance_transform(feature_transform(binarized))
@. dist = 1 - dist

segments = watershed(dist, dist .< distance_threshold)
#dist =  distance_transform(feature_transform(rescaled.im .> level_threshold))
#markers = label_components(dist .< distance_threshold)
return map(i->get_random_color(i), labels_map(segments)) .* (1 .-binarized)  #labels_map(segments) #.* (1 .-binarized)
end

using .ThermovisorImages

im = read_temperature_file(raw"D:\JuliaDepoth\dev\ThermovisorData.jl\thermal images\T500.csv")
roi = ThermovisorImages.fit_all_patterns(im)


    using ImageShow
    source_dir = joinpath(abspath(joinpath(@__DIR__, "..")),"src")

    include(joinpath(source_dir,"ThermovisorImages.jl"))
    im1 = ThermovisorImages.read_temperature_file(raw"E:\JULIA\JULIA_DEPOT\dev\ThermovisorData\thermal images\T500.csv")
    simshow(im1.initial)
    #c = ThermovisorImages.CircleObj()
    #fltrd = ThermovisorImages.filter_image(im1;label = 1)
    #ThermovisorImages.fit_centred_obj!(c,fltrd)
    #markers = ThermovisorImages.marker_image(im1)
    roi = ThermovisorImages.fit_all_patterns(im1,ThermovisorImages.CircleObj)
    if isempty(roi)
        return "empty"
    end
    y = copy(im1.initial)
    y_view = ThermovisorImages.c_view(y,roi[1])
    ThermovisorImages.recalculate_with_new_emissivity!(y_view,0.5,1.0,λ_left=14.0, λ_right=nothing) #single wavelength version
    simshow(y)
    maximum(y)
    y2 = copy(im1.initial)
    y2_view = ThermovisorImages.c_view(y2,roi[1])
    ThermovisorImages.recalculate_with_new_emissivity!(y2,roi[1],0.5,1.0,λ_left=14.0, λ_right=nothing)
    maximum(im1.initial)
    maximum(y2)

    using BenchmarkTools

    @benchmark ThermovisorImages.recalculate_with_new_emissivity!(y_view,0.5,1.0,λ_left=14.0, λ_right=nothing)