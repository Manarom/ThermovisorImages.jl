using Revise,ImageShow,Images,BenchmarkTools
source_dir = joinpath(abspath(joinpath(@__DIR__, "..")),"src")
includet(joinpath(source_dir,"ThermovisorData.jl"))
begin 
img = load(download("http://docs.opencv.org/3.1.0/water_coins.jpg"));

im_float = Float64.(Gray.(img))
bw = im_float .> 0.5;

dist = 1 .- distance_transform(feature_transform(bw));

markers = label_components(dist .< -15);

segments = watershed(dist, markers)
simshow(labels_map(segments))
markers2 = labels_map(segments) .* (1 .-bw) 
simshow(markers2 )  

im_rescaled = ThermovisorData.RescaledImage(1 .-im_float)
markers = ThermovisorData.marker_image(im_rescaled)
simshow(markers.markers)
end
begin 
    (im_rescaled2,) = ThermovisorData.read_temperature_file(raw".\thermal images\T600.csv")

    markers2 = ThermovisorData.marker_image(im_rescaled2)
    simshow(markers2.markers)

    #maximum(markers2)
end
typeof(markers2)
extflag = ThermovisorData.external_flag(markers2,1)
simshow(extflag)
im_copy = copy(im_rescaled2.im)
v=  @view im_copy[extflag]
typeof(v)
ThermovisorData.filter_image!( im_copy,extflag)