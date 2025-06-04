### A Pluto.jl notebook ###
# v0.20.6

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ fc6af4b0-1127-11f0-1b66-a59d87c5b141
begin# we need the project structrue inside the notebook
	notebook_path = @__DIR__()# this notebook local path
	project_path = abspath(joinpath(notebook_path,".."))#project's path
	import Pkg
	Pkg.activate(project_path)
	sources_path = joinpath(project_path,"src")# it is supposed that sources are in separate folder \project_folder\src\
	assets_folder = joinpath(project_path,"assets")
end;

# ╔═╡ 051044c5-760c-4b60-90fc-82a347c3b6bc
using Revise,PlutoUI,LaTeXStrings,Images,ImageShow,Plots,BenchmarkTools,Dates,FileIO,ImageIO,Optim

# ╔═╡ e71f123c-284e-4107-8231-4031873f122c
using Main.ThermovisorData # it always fails on the first run when Revise is in  use 

# ╔═╡ 4460f260-f65f-446d-802c-f2197f4d6b27
md"""
### This notebook demostrates main features of `ThermovisorData.jl` package
#
**ThervisorData.jl** is a small package designed to process static thermal images stored as matrices. Each matrix element represents a temperature value. This package enables users to load images from csv files, calculate temperature distributions, and compute statistical analyses for temperatures along specified lines. It also calculates averaged angular and radial temperature distributions (along with standard deviations) within Regions of Interest (ROIs) such as circles, squares, and rectangles. These ROI objects can be fitted to thermally distinct areas (regions that stand out from the background).

The image can also be filtered to remove all other patterns from the image except the one with selected mark.

This package was written to study the temperature distribution across the heated sample for the emissivity measuring facility described in this [`paper`](https://link.springer.com/article/10.1007/s00340-024-08331-9)

"""

# ╔═╡ 2c5e6e4c-92af-4991-842a-7e5bdc55a46d
md"""
###### This notebooks describes the following features:

*  1. Loading image
*  2. Creating `CentredObj` marker (ROI object)
*  3. Filtering the image
*  4. Fitting marker position and size to the pattern
*  5. Evaluating radial and angular temperature distributions
*  6. Fitting multiple ROI objects to the image with several temperature features
"""

# ╔═╡ 215ed2f4-71ba-4cb5-b198-677d0d7ffb38
md" default image saving folder $(@bind image_save_folder PlutoUI.TextField(default = assets_folder))"

# ╔═╡ f6c1be87-94d2-4b08-a52d-6eb637192ee8
includet(joinpath(sources_path,"ThermovisorData.jl"))
# include(joinpath(sources_path,"ThermovisorData.jl")) #replace includet with include if Revise is not needed

# ╔═╡ 870113c3-b439-4d34-90d8-fdd8a158f9dd
md"""
#### 1. Loading image


All examples of temperature distribution images from real application are in 
`"...project_folder\thermal images\ "` folder. All thermal images are in csv - format

Function `ThermovisorData.find_temperature_files(images_folder)` returns `Dict` with thermovisor data file labels matched to the full file names and actual measurements temperatures. This works simple by parsing the name of the file, to be parsed as thermal image file it's name should be like :
"any-name-T567.csv" . This file is automatically interpreted as a thermal image recorded for the temperature T = 567, here all numbers after the "T" symbol are interpreted as the temperature value.

"""

# ╔═╡ c241bae6-1925-4be1-af41-44673f02617a
begin 
	images_folder = joinpath(project_path,"thermal images")
	files_in_dir = ThermovisorData.find_temperature_files(images_folder)
end

# ╔═╡ 43a1fb58-cd5e-4634-8770-0ff1809b2191
md"Select file tag: $(@bind temp_to_load Select([v for v in keys(files_in_dir)])) "

# ╔═╡ 794ebd5e-e9e0-4772-98a9-43e20c7ef4da
#reading selected file
begin 
	(rescaled_image,file_date) = read_temperature_file(files_in_dir[temp_to_load][2])
	Dates.unix2datetime(file_date)
end

# ╔═╡ 429cf33f-4422-44f0-beb8-5a1908a72273
md"""
	### 2.Creating `CentredObj` marker

	`CentredObj` acts as a region of interest (ROI), it has independent (from the image) coordinates and size, thus it can be used to scan the image and monitor it's values. `CentredObj` ROI can be of different shapes viz circle, rectangle and square. 

	Select marker type and adjust it's centre coordinates and dimentions, use sliders to move the ROI over the image:	

	marker type = $(@bind test_mask_type Select([CircleObj=>"circle",  SquareObj =>"square",RectangleObj =>"rectangle"])) \
	center point row index = $(@bind test_mask_center_x Slider(1:1000,default=200,show_value=true)) \
	center point column index = $(@bind test_mask_center_y Slider(1:1000,default=200,show_value=true)) \
	side a= $(@bind test_side_a Slider(1:1000,default=150,show_value=true)) 
	side b= $(@bind test_side_b Slider(1:1000,default=100,show_value=true))

	show marker on separate image = $(@bind is_make_obj_selfy CheckBox(false))

	"""

# ╔═╡ 854731c1-7a34-4066-aa74-01629c87d75d
begin

	mm_per_pixel = Ref(0.8); # this is used to calibrate the image dimentions 
	test_coord_vect = [test_mask_center_x,test_mask_center_y,test_side_a]
	if test_mask_type==RectangleObj
		append!(test_coord_vect,test_side_b)
	end
	test_centre_obj = ThermovisorData.obj_from_vect(test_mask_type,test_coord_vect)
	is_make_obj_selfy ? centre_obj_image = draw(test_centre_obj,thickness=25) : nothing
	
	rgb_image_initial = draw!(rescaled_image.initial,color_scheme="HEAT",test_centre_obj,show_cross = true)
			imag_initial_show = ThermovisorData.draw_line_within_mask!(rgb_image_initial,test_centre_obj,0,10,color_scheme="R1")
	
end

# ╔═╡ 48c53ed9-127d-4e8e-bd29-416292057bff
ThermovisorData.convert_to_drawable(test_centre_obj)

# ╔═╡ 9f55febb-b047-4b22-8575-209d45354d51
md" Export image to file $(@bind save_image CheckBox(default = false))"

# ╔═╡ 4feea216-ee48-42a3-b4ba-454f28ff690a
begin
	if save_image 		
		FileIO.save(joinpath(image_save_folder,"initial_image.png"), rgb_image_initial)
	end
end

# ╔═╡ 13f01881-2645-429b-9856-6c3f19c0ad48
md"""
The following block evaluates the average and standart deviation of temperature within the ROI 

Average temperature within the marker <T>=$(ThermovisorData.mean_within_mask(rescaled_image.initial,test_centre_obj)) 

Standat deviation of average std(T) = $(ThermovisorData.std_within_mask(rescaled_image.initial,test_centre_obj))
"""

# ╔═╡ 5a212007-c0e8-4b1b-94d1-30bdb1efdb9c
md"""
### 2. Filtering image
There are several methods for `ThermovisorData.filter_image` function. When the image is the only input argument, it goes through several operation from `ImageSegmentation` package.

*  1.Binarization
*  2.Distance transform
*  3.Labeling 
*  4.Filtering by selecting one of the labels

After labeling, by default it takes the last label (maximum label value), but label value could be provided externally with a corresponding keyword `label`
If the `CnetredObj` is provided as a secon input argument filtering just removes all elements of the initial image which are not within the `CentredObj`.

"""

# ╔═╡ 3fbc6b45-974e-430e-a4e6-960323015e74
md""" 

do you want to filter image ? = $(@bind filter_init_image CheckBox(default=true))

filter by $(@bind filter_by_option Select(["CentredObj", "Pattern"]))

show reduced $(@bind is_show_filtered_reduced CheckBox(default=false))

"""

# ╔═╡ ca5eea20-2bb3-4407-aa09-af8de2332b84
md"##### Filtered image:"

# ╔═╡ 8b6f604d-157b-42cd-a0c6-8bd5562b47ef
begin 
	if filter_init_image
		# filtering the image
		if filter_by_option=="CentredObj"
			filtered_by_obj = ThermovisorData.filter_image(rescaled_image,test_centre_obj)
		else
			filtered_by_obj = ThermovisorData.filter_image(rescaled_image)
		end
		h2 = ThermovisorData.draw(filtered_by_obj,draw_reduced=is_show_filtered_reduced) 
	end
end

# ╔═╡ 38a45961-0ffb-43d4-aa24-36d503ed4618
md"Save the current heatmap $(@bind save_distr CheckBox(default=false))"

# ╔═╡ 1467b184-22ac-4038-ad1b-f084d4443b27
save_distr ? savefig(joinpath(notebook_path,"heatmap.png")) : nothing

# ╔═╡ c87a830a-f48a-4444-81bc-3efd69a130ad
md"""

###  4,5. Fitting marker and Evaluating radial and angular temperature distributions

This segment illustrates an example of temperature distribution within the Region of Interest (ROI) calculation, as well as how the heated object is fitted within this analysis.

"""

# ╔═╡ 6d37916c-7895-49d3-b8a3-c8661050ebcb
md"""
There are two algorithms available for drawing a line within an image. One of these is the `xiaolin_wu` function from the `ImageDraw` package. This algorithm produces two coordinates for every point along the line, and the resulting temperature for each position is calculated as the average of these two points. Another algorithm is based on `ImageDraw.bresenham` it is used by default.

Try Xiaolin-Wu algorithm? $(@bind is_use_wu CheckBox(default=false))

"""

# ╔═╡ 6482d05d-06e2-43cc-ab53-ff4bbcd63e3e
md"mm per pixel calibration value = $(mm_per_pixel[])"

# ╔═╡ c67290fc-6291-4f3e-a660-a3c4afa3a5e3
md"""
	Creating mask object:
	image type = $(@bind image_type Select(["filtered", "filtered full","not filtered"])) \
	ROI type = $(@bind mask_type Select([CircleObj=>"circle",  SquareObj =>"square",RectangleObj =>"rectangle"])) \
	center point x= $(@bind mask_center_x Slider(1:1000,default=200,show_value=true)) \
	center point y= $(@bind mask_center_y Slider(1:1000,default=200,show_value=true)) \
	side a= $(@bind side_a Slider(1:1000,default=150,show_value=true)) \
	side b= $(@bind side_b Slider(1:1000,default=150,show_value=true))

	"""

# ╔═╡ 4e1a5050-59b0-4d24-98bb-1520c06b28c5
begin 
	coord_vect = [mask_center_y,mask_center_x,side_a]
	if mask_type==RectangleObj
		append!(coord_vect,side_b)
	end
	centre_obj = ThermovisorData.obj_from_vect(mask_type,coord_vect)
	#is_make_obj_selfy ? centre_obj_image = draw(centre_obj,thickness=25) : nothing
end

# ╔═╡ 71eb240a-5a45-4bf3-b35c-a5820ca6da6c
md" Hitting this checkbox forces the ROI to fit the temperature pattern $(@bind is_fit_mask CheckBox(false))"

# ╔═╡ 768535e0-a514-4dff-ac8b-0d7ca126149c
# fitting the loaded image
begin 

	filtered_by_markers = filter_image(rescaled_image,marker_image(rescaled_image))
	fitted_obj = ThermovisorData.copyobj(centre_obj)	
	if image_type=="filtered"
		image_to_show = ThermovisorData.reduced_image(filtered_by_markers)
	elseif image_type =="not filtered"
		image_to_show = rescaled_image.initial
	else # full filtered
		image_to_show = filtered_by_markers.full.initial
	end
	if is_fit_mask 
		fit_centred_obj!(fitted_obj,filtered_by_markers,fit_reduced = image_type=="filtered") 
		@show fitted_obj
		mm_per_pixel[]=sqrt((0.25π*25^2)/area(fitted_obj))
	else
		mm_per_pixel[]=1.0
	end
	
end;

# ╔═╡ 5d6222cf-99f3-4ce9-a4a2-91c17dc9c0d2
fitted_obj

# ╔═╡ 59f9a7f2-9601-431c-a897-543fa25c64c4
fitted_obj

# ╔═╡ e1ccfd33-3d54-4249-86f1-381a1ef90615
md"""
The upper figure shows the ROI and the inclined line which goes through its center. By adjusting ROI position, the orientation and the length of this line temperature distribution of some feature can be studied. \
 	direction angle in degrees = $(@bind direction_angle Slider(0:1:360,show_value=true,default=45)) \
	line length in $(is_fit_mask ? "mm" : "pixels" )= $(@bind line_length Slider(0.1:0.1:250,show_value=true,default=100))
"""

# ╔═╡ 42a7b186-aa04-4249-a129-bf925f181008
begin
	rgb_image = draw!(image_to_show,color_scheme="HEAT",fitted_obj,show_cross = true)
	imag = ThermovisorData.draw_line_within_mask!(rgb_image,fitted_obj,direction_angle,line_length/mm_per_pixel[],color_scheme="R1")
end

# ╔═╡ b096c4f2-9dce-409d-874a-a851f577bf92
begin 
	#ThermovisorData.diag_ang(fitted_obj)
	(along_line_length,distrib,line_points) = ThermovisorData.along_mask_line_distribution(image_to_show,fitted_obj,direction_angle,line_length,use_wu=is_use_wu,length_per_pixel=mm_per_pixel[])
	#@show length(points)
	
	pl_distrib=ThermovisorData.plot_along_line_distribution(along_line_length,distrib,is_centered = true)
	#pl_distrib
end

# ╔═╡ d45f106d-032a-4102-a26f-7393c2220f72
md"Evaluate radial distribution $(@bind is_rad_distr_eval CheckBox(false))"

# ╔═╡ 39e50296-21ff-4407-894f-2a380dc51e21
begin 
	
	if is_rad_distr_eval
		new_line =line_length# minimum(side(fitted_obj))-5
		ang_range = 0.0:1:60 # range of angles 
		R,D = ThermovisorData.radial_distribution(image_to_show,fitted_obj,ang_range,line_length=new_line,length_per_pixel=mm_per_pixel[])

		(L,meanD,stdD,l_b,u_b,t_values) = ThermovisorData.radial_distribution_statistics(R,D)

		
		p_radial = ThermovisorData.plot_radial_distribution_statistics(L,meanD,stdD,
        	l_b,u_b)

	end
end

# ╔═╡ ca05bd4f-5656-4531-b357-331c62661174
begin 
	if is_rad_distr_eval
	angs = collect(ang_range)
	(angs_plot,Dang,stdDang,l_bang,u_bang,t_valuesang) = ThermovisorData.angular_distribution_statistics(angs,R,D)
	
	p_angular = ThermovisorData.plot_angular_distribution_statistics(angs_plot,Dang,stdDang,
        	l_bang,u_bang)
	end
end

# ╔═╡ e9216d7a-c2f3-44c0-a7d9-2c62ac35ecd9
md"Save image with marker and temperature distributions $(@bind save_average_radial_distribution CheckBox(default=false))"


# ╔═╡ b4ce12e3-29ec-41ac-89d3-06d08ef2beca
begin 
	if is_rad_distr_eval&& save_average_radial_distribution  	
		FileIO.save(joinpath(assets_folder,"filtered_image_with_marker.png"),imag)
		savefig(pl_distrib,joinpath(assets_folder,"line_distrib.png"))
		savefig(p_radial,joinpath(assets_folder,"radial_distrib.png"))
		savefig(p_angular,joinpath(assets_folder,"angular_distrib.png"))
	end
end;

# ╔═╡ 764a320c-ff6b-48d0-a5b4-48a3df3ece01
md"""
###  6. Fitting multiple ROI objects to the image with several temperature features

In this block we are going to create the image with multiple randomly distributed  patterns and fit a vector of  `CentredObj` ROIs to this image by calling `ThermovisorData.fit_all_patterns` function.

"""

# ╔═╡ 10954f10-9414-4839-872f-c2516d5d8e4e
md"""
	Click this button to generate the image pattern $(@bind fit_multiple Button("Regenerate patterns"))
	"""

# ╔═╡ 6e728ea6-38be-437a-96b4-9fa084f8fec5
md"""
Select objects ROI type = $(@bind multifit_roi_type Select([CircleObj=>"circle",  SquareObj =>"square",RectangleObj =>"rectangle"])) 
"""

# ╔═╡ cc909b53-ed4d-44a1-a410-ff25533afc2d
md"Initial image with randomly distributed patterns"

# ╔═╡ d5b6f453-5e92-41e6-a45f-cb75660bc198
begin # generating pattern 
	fit_multiple
	Patterns_number = 10 #total number of patterns
	image_size = (480,640)
	img = fill(0.0,image_size...)# filling initial scene
	rnd_centre() = [rand(1:image_size[1]),rand(1:image_size[2])] # random center positions
	rnd_diam() = rand(10:10:30) # random diameter generator
	# filling the vector of initial ROIs
	centered_objs = [CircleObj(rnd_centre(),rnd_diam()) for _ in 1:Patterns_number]
	for c in centered_objs
		img[c]=1 # centred objects can be used as indices in images 
	end
	rs = RescaledImage(img)
	markers = ThermovisorData.marker_image(rs,distance_threshold = 0.0)
	separate_patterns_number = length(markers) # now we need to check for the separate patterns number 
	rgb_markers = ThermovisorData.draw(Float64.(markers.markers))	# converting to rgb image
end

# ╔═╡ 96ad6e27-52dd-41aa-b115-f852049a485a
md"""Number of separate patterns:    $(separate_patterns_number)"""

# ╔═╡ 0badf26a-38fa-45be-9704-d4e80b12a9cb
md"""
	Select this checkbox to fit all objects $(@bind is_fit_multiple CheckBox(default=false))
	"""

# ╔═╡ 6adfae4d-5137-4692-b9f3-3793c4c76202
begin # fitting ROI's to image with several 
	fit_multiple
	if is_fit_multiple

		fitted_rois = ThermovisorData.fit_all_patterns(rs,multifit_roi_type,distance_threshold = 0.0,sort_by_area=true)
		rgb_image_multi_roi = draw!(img,fitted_rois[1],show_cross = true,fill=true)
		for i in 2:length(fitted_rois)
			 ThermovisorData.draw!(rgb_image_multi_roi,fitted_rois[i],thickness = 1,fill=true,show_cross = true)
		end
		rgb_image_multi_roi
		#
	else
		nothing
	end
end

# ╔═╡ 68b33b39-5ef5-4560-b4b2-1fe2f43a3628
md" Save multiple patterns fit ? $(@bind is_save_multipattern_fit CheckBox(default=false))"

# ╔═╡ 8a132aba-aa8a-428a-84a2-0ab6e5e2b891
begin 
	if is_save_multipattern_fit
		FileIO.save(joinpath(assets_folder,"multiple_patterns_initial.png"),rgb_markers)
		if is_fit_multiple
			FileIO.save(joinpath(assets_folder,"multiple_patterns_fitted.png"),rgb_image_multi_roi)
		end
	end
end

# ╔═╡ 821a7c95-f4da-410d-b780-111abb6d0db5
md"""
	Show/hide rois $(@bind is_draw_rois Select( ["show" ;"hide" ]))
	"""

# ╔═╡ 2af862c2-c026-43a8-82e3-77ff8cb2095c
begin
	im_coin = load(download("http://docs.opencv.org/3.1.0/water_coins.jpg"))
	im_coin_float = 1 .- Float64.(Gray.(im_coin))
	
	rescaled_coin = ThermovisorData.RescaledImage(im_coin_float)
	markers_coins = ThermovisorData.marker_image(rescaled_coin)
	
	fitted_rois_coins = ThermovisorData.fit_all_patterns(rescaled_coin,multifit_roi_type)
	
	rgb_image_coins = draw!(im_coin_float,fitted_rois_coins[1],show_cross = true,fill=true)
	for i in 2:length(fitted_rois_coins)
			 ThermovisorData.draw!(rgb_image_coins,fitted_rois_coins[i],thickness = 1,fill=true,show_cross = true)
	end
	rgb_image_coins
end

# ╔═╡ b640fcd0-3e49-471d-b281-87137a781eba
begin 
	test_markers = ThermovisorData.marker_image(rescaled_coin)
	before_sorting = simshow(test_markers.markers)
	inds = Vector{Int}(undef,length(test_markers))
	@show areas_vals = ThermovisorData.areas(test_markers)
	sortperm!(inds,areas_vals,rev=true)
	permute!(test_markers.ViewsVect,inds)
	@show ThermovisorData.areas(test_markers)
	for i = 1:length(test_markers)
        #v = view(test_markers.markers,test_markers.ViewsVect[i])
		#fill!(v,i)
		test_markers[i] = i
    end
	@show ThermovisorData.areas(test_markers)
	after_sorting = simshow(test_markers.markers)
end;

# ╔═╡ 7a00ce43-94e2-4f68-b651-b57bf7d6ab05
hcat(before_sorting,after_sorting)

# ╔═╡ ecfd5449-29f6-452b-ae3d-d8e40932c8a0
simshow(ThermovisorData.shrinked_flag(test_markers,1))

# ╔═╡ 1fb56e7a-6bc5-4c9d-9de1-22a62e39bb66
c1 = ThermovisorData.CircleObj()

# ╔═╡ a7d6f66c-0774-419b-85ee-b8da74498227
fl1 = ThermovisorData.shrinked_flag(test_markers,1);

# ╔═╡ 3e408b86-bb21-4d48-8842-a21b054726be
x0 = [1.0; 2.0; 3.0]

# ╔═╡ 9d4a3283-26f4-45fb-ab26-7f55d9f01a4e
ThermovisorData.fill_x0!(x0,fl1,c1)

# ╔═╡ 4adf3ee6-742d-421d-bb54-dca02aebd833
ThermovisorData.fit_centred_obj!(c1,fl1)

# ╔═╡ 81b1776f-c527-4235-b84e-c9c912098220


# ╔═╡ 3ab979f6-3a5f-4739-9c49-ff1e75a94a62
ThermovisorData.area(c1)

# ╔═╡ 0903ce52-3fe8-41fa-a0d8-bc65fa589d10
begin 
	#cent23 = RectangleObj([48,50],[35,25])
	#im23= ThermovisorData.cent_to_flag(cent23,(100,100))
	im23 = image_to_show .>30
	fff = findall(im23)
	 (min_i,max_i) = extrema(fff)
	dinds = max_i - min_i + CartesianIndex(1,1)#indices difference
	
	fff_reduced = fill(false,Tuple.(dinds))
	for inds in fff
		#@show inds
		inds_in = inds - min_i + CartesianIndex(1,1)
		#@show inds_in
		fff_reduced[inds_in] = im23[inds]
	end
	
end

# ╔═╡ Cell order:
# ╟─4460f260-f65f-446d-802c-f2197f4d6b27
# ╟─2c5e6e4c-92af-4991-842a-7e5bdc55a46d
# ╠═fc6af4b0-1127-11f0-1b66-a59d87c5b141
# ╠═051044c5-760c-4b60-90fc-82a347c3b6bc
# ╟─215ed2f4-71ba-4cb5-b198-677d0d7ffb38
# ╠═f6c1be87-94d2-4b08-a52d-6eb637192ee8
# ╠═e71f123c-284e-4107-8231-4031873f122c
# ╟─870113c3-b439-4d34-90d8-fdd8a158f9dd
# ╠═c241bae6-1925-4be1-af41-44673f02617a
# ╟─43a1fb58-cd5e-4634-8770-0ff1809b2191
# ╟─794ebd5e-e9e0-4772-98a9-43e20c7ef4da
# ╟─429cf33f-4422-44f0-beb8-5a1908a72273
# ╟─854731c1-7a34-4066-aa74-01629c87d75d
# ╟─48c53ed9-127d-4e8e-bd29-416292057bff
# ╟─9f55febb-b047-4b22-8575-209d45354d51
# ╟─4feea216-ee48-42a3-b4ba-454f28ff690a
# ╟─13f01881-2645-429b-9856-6c3f19c0ad48
# ╟─5a212007-c0e8-4b1b-94d1-30bdb1efdb9c
# ╟─3fbc6b45-974e-430e-a4e6-960323015e74
# ╟─ca5eea20-2bb3-4407-aa09-af8de2332b84
# ╟─8b6f604d-157b-42cd-a0c6-8bd5562b47ef
# ╟─38a45961-0ffb-43d4-aa24-36d503ed4618
# ╟─1467b184-22ac-4038-ad1b-f084d4443b27
# ╟─c87a830a-f48a-4444-81bc-3efd69a130ad
# ╟─768535e0-a514-4dff-ac8b-0d7ca126149c
# ╟─5d6222cf-99f3-4ce9-a4a2-91c17dc9c0d2
# ╟─6d37916c-7895-49d3-b8a3-c8661050ebcb
# ╟─6482d05d-06e2-43cc-ab53-ff4bbcd63e3e
# ╟─c67290fc-6291-4f3e-a660-a3c4afa3a5e3
# ╟─4e1a5050-59b0-4d24-98bb-1520c06b28c5
# ╠═42a7b186-aa04-4249-a129-bf925f181008
# ╟─e1ccfd33-3d54-4249-86f1-381a1ef90615
# ╟─b096c4f2-9dce-409d-874a-a851f577bf92
# ╟─59f9a7f2-9601-431c-a897-543fa25c64c4
# ╟─39e50296-21ff-4407-894f-2a380dc51e21
# ╟─ca05bd4f-5656-4531-b357-331c62661174
# ╟─71eb240a-5a45-4bf3-b35c-a5820ca6da6c
# ╟─d45f106d-032a-4102-a26f-7393c2220f72
# ╟─e9216d7a-c2f3-44c0-a7d9-2c62ac35ecd9
# ╟─b4ce12e3-29ec-41ac-89d3-06d08ef2beca
# ╟─764a320c-ff6b-48d0-a5b4-48a3df3ece01
# ╠═10954f10-9414-4839-872f-c2516d5d8e4e
# ╟─6e728ea6-38be-437a-96b4-9fa084f8fec5
# ╠═cc909b53-ed4d-44a1-a410-ff25533afc2d
# ╟─d5b6f453-5e92-41e6-a45f-cb75660bc198
# ╟─96ad6e27-52dd-41aa-b115-f852049a485a
# ╟─0badf26a-38fa-45be-9704-d4e80b12a9cb
# ╠═6adfae4d-5137-4692-b9f3-3793c4c76202
# ╟─68b33b39-5ef5-4560-b4b2-1fe2f43a3628
# ╟─8a132aba-aa8a-428a-84a2-0ab6e5e2b891
# ╠═821a7c95-f4da-410d-b780-111abb6d0db5
# ╠═2af862c2-c026-43a8-82e3-77ff8cb2095c
# ╠═b640fcd0-3e49-471d-b281-87137a781eba
# ╠═7a00ce43-94e2-4f68-b651-b57bf7d6ab05
# ╠═ecfd5449-29f6-452b-ae3d-d8e40932c8a0
# ╠═1fb56e7a-6bc5-4c9d-9de1-22a62e39bb66
# ╠═a7d6f66c-0774-419b-85ee-b8da74498227
# ╠═3e408b86-bb21-4d48-8842-a21b054726be
# ╠═9d4a3283-26f4-45fb-ab26-7f55d9f01a4e
# ╠═4adf3ee6-742d-421d-bb54-dca02aebd833
# ╠═81b1776f-c527-4235-b84e-c9c912098220
# ╠═3ab979f6-3a5f-4739-9c49-ff1e75a94a62
# ╠═0903ce52-3fe8-41fa-a0d8-bc65fa589d10
