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
end;

# ╔═╡ 051044c5-760c-4b60-90fc-82a347c3b6bc
using Revise,PlutoUI,LaTeXStrings,Images,ImageShow,Plots,BenchmarkTools,Dates,Images

# ╔═╡ e71f123c-284e-4107-8231-4031873f122c
using Main.ThermovisorData

# ╔═╡ 4460f260-f65f-446d-802c-f2197f4d6b27
md"""
### This notebook demostrates main features of ThermovisorData package
#
**ThervisorData.jl** is a small package designed to process static thermal images stored as matrices in CSV format, where each matrix element represents a temperature value. This package enables users to load images from files, calculate temperature distributions, and compute statistical analyses for temperatures along specified lines. It also calculates averaged angular and radial temperature distributions (along with standard deviations) within Regions of Interest (ROIs) such as circles, squares, and rectangles. These ROI objects can be fitted to thermally distinct areas (relative to their surroundings), such as the most heated regions within the scene.

The image can also be filtered to remove all other patterns from the image except the one with selected mark.

This package was written to study the temperature distribution across the heated sample for the emissivity measuring facility described in this [`paper`](https://link.springer.com/article/10.1007/s00340-024-08331-9)

"""

# ╔═╡ 2c5e6e4c-92af-4991-842a-7e5bdc55a46d
md"""
###### This notebooks describes the following features:

*  1. Loading image
*  2. Creating `CentredObj` marker
*  3. Filtering the image
*  4. Fitting marker position and size to the pattern
*  5. Evaluating radial and angular temperature distributions
*  6. Fitting multiple ROI objects to the image with several temperature features
"""

# ╔═╡ f6c1be87-94d2-4b08-a52d-6eb637192ee8
includet(joinpath(sources_path,"ThermovisorData.jl"))

# ╔═╡ 870113c3-b439-4d34-90d8-fdd8a158f9dd
md"""
#### 1. Loading image


All examples of temperature distribution images from real application are in 
`"project_folder\thermal images\ "` folder. All thermal images are in csv - format

Function `ThermovisorData.find_temperature_files(images_folder)` returns `Dict` with tthermovisor data file labels matched to the full file names and actual measurements temperatures. This works simple by parsing the name of the file, it should contain :
"any-name-T567.csv" is interpreted as a thermal image recorded for the temperature T = 567, here all numbers after the "T" symbol are interpreted as the temperature value.

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

	Select marker type and adjust it's centre coordinates and dimentions
	
	marker type = $(@bind test_mask_type Select([CircleObj=>"circle",  SquareObj =>"square",RectangleObj =>"rectangle"])) \
	center point row index = $(@bind test_mask_center_x Slider(1:1000,default=100,show_value=true)) \
	center point column index = $(@bind test_mask_center_y Slider(1:1000,default=100,show_value=true)) \
	side a= $(@bind test_side_a Slider(1:1000,default=100,show_value=true)) 
	side b= $(@bind test_side_b Slider(1:1000,default=100,show_value=true))

	automake image = $(@bind is_make_obj_selfy CheckBox(false))

	"""

# ╔═╡ a963a19c-1a1f-4f44-975b-7803a6b9d9cd
begin 
	test_coord_vect = [test_mask_center_x,test_mask_center_y,test_side_a]
	if test_mask_type==RectangleObj
		append!(test_coord_vect,test_side_b)
	end
	test_centre_obj = ThermovisorData.obj_from_vect(test_mask_type,test_coord_vect)
	is_make_obj_selfy ? centre_obj_image = draw(test_centre_obj,thickness=25) : nothing
end

# ╔═╡ 9ffe9b99-ab18-49b2-9179-3afb81b7be48
mm_per_pixel = Ref(0.8);

# ╔═╡ 5a212007-c0e8-4b1b-94d1-30bdb1efdb9c
md"""
### 2. Filtering image
*  2.Binarization
*  3.Distance transform
*  4.Labeling image
*  5.Filtering image with one label
"""

# ╔═╡ 3fbc6b45-974e-430e-a4e6-960323015e74
md""" do you want to filter initial image ? = $(@bind filter_init_image CheckBox(default=false))"""

# ╔═╡ 854731c1-7a34-4066-aa74-01629c87d75d
begin
	if filter_init_image
		filtered_by_obj = filter_image(rescaled_image,test_centre_obj)
		h2 = heatmap(filtered_by_obj.full.initial[end:-1:1,:],dpi=600)
		title!("\\mu =$(filtered_mean(filtered_by_obj)); std = $(filtered_std(filtered_by_obj))")
	else
			rgb_image_initial = draw!(rescaled_image.initial,color_scheme="HEAT",test_centre_obj,show_cross = true)
			imag_initial_show = ThermovisorData.draw_line_within_mask!(rgb_image_initial,test_centre_obj,0,10,color_scheme="R1")
	end
end

# ╔═╡ 13f01881-2645-429b-9856-6c3f19c0ad48
md"""
Average temperature within the marker <T>=$(ThermovisorData.mean_within_mask(rescaled_image.initial,test_centre_obj)) 

Standat deviation of average std(T) = $(ThermovisorData.std_within_mask(rescaled_image.initial,test_centre_obj))
"""

# ╔═╡ 38a45961-0ffb-43d4-aa24-36d503ed4618
md"Save the current heatmap $(@bind save_distr CheckBox(default=false))"

# ╔═╡ 1467b184-22ac-4038-ad1b-f084d4443b27
save_distr ? savefig(joinpath(notebook_path,"heatmap.png")) : nothing

# ╔═╡ c87a830a-f48a-4444-81bc-3efd69a130ad
md"""

	### Testing work with in-marker distribution
	"""

# ╔═╡ 6d37916c-7895-49d3-b8a3-c8661050ebcb
md" Use Wu algorithm? $(@bind is_use_wu CheckBox(default=false))"

# ╔═╡ 6482d05d-06e2-43cc-ab53-ff4bbcd63e3e
md"mm per pixel calibration value = $(mm_per_pixel[])"

# ╔═╡ c67290fc-6291-4f3e-a660-a3c4afa3a5e3
md"""
	Creating mask object:

	image type = $(@bind image_type Select(["filtered", "filtered full","not filtered"]))
	mask type = $(@bind mask_type Select([CircleObj=>"circle",  SquareObj =>"square",RectangleObj =>"rectangle"])) \
	center point x= $(@bind mask_center_x Slider(1:1000,default=100,show_value=true)) \
	center point y= $(@bind mask_center_y Slider(1:1000,default=100,show_value=true)) \
	side a= $(@bind side_a Slider(1:1000,default=100,show_value=true)) \
	side b= $(@bind side_b Slider(1:1000,default=100,show_value=true))

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
md" Fit the mask ? $(@bind is_fit_mask CheckBox(false))"

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

# ╔═╡ e1ccfd33-3d54-4249-86f1-381a1ef90615
md"""
 	direction angle in degrees = $(@bind direction_angle Slider(0:1:360,show_value=true,default=45))
	line length in $(is_fit_mask ? "mm" : "pixels" )= $(@bind line_length Slider(0.1:0.1:250,show_value=true,default=20))
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
		ang_range = 0.0:1:180
		R,D = ThermovisorData.radial_distribution(image_to_show,fitted_obj,ang_range,line_length=new_line,length_per_pixel=mm_per_pixel[])

		(L,meanD,stdD,l_b,u_b,t_values) = ThermovisorData.radial_distribution_statistics(R,D)

		
		ThermovisorData.plot_radial_distribution_statistics(L,meanD,stdD,
        	l_b,u_b)

	end
end

# ╔═╡ ca05bd4f-5656-4531-b357-331c62661174
begin 
	if is_rad_distr_eval
	angs = collect(ang_range)
	(angs_plot,Dang,stdDang,l_bang,u_bang,t_valuesang) = ThermovisorData.angular_distribution_statistics(angs,R,D)
	
	ThermovisorData.plot_angular_distribution_statistics(angs_plot,Dang,stdDang,
        	l_bang,u_bang)
	end
end

# ╔═╡ e9216d7a-c2f3-44c0-a7d9-2c62ac35ecd9
md"Save averaged radial temperature ditribution $(@bind save_average_radial_distribution CheckBox(default=false))"


# ╔═╡ b4ce12e3-29ec-41ac-89d3-06d08ef2beca
save_average_radial_distribution ? savefig(p,joinpath(notebook_path,"radial_distrib.png")) : nothing

# ╔═╡ Cell order:
# ╟─4460f260-f65f-446d-802c-f2197f4d6b27
# ╠═2c5e6e4c-92af-4991-842a-7e5bdc55a46d
# ╠═fc6af4b0-1127-11f0-1b66-a59d87c5b141
# ╠═051044c5-760c-4b60-90fc-82a347c3b6bc
# ╠═f6c1be87-94d2-4b08-a52d-6eb637192ee8
# ╠═e71f123c-284e-4107-8231-4031873f122c
# ╠═870113c3-b439-4d34-90d8-fdd8a158f9dd
# ╟─c241bae6-1925-4be1-af41-44673f02617a
# ╠═43a1fb58-cd5e-4634-8770-0ff1809b2191
# ╟─794ebd5e-e9e0-4772-98a9-43e20c7ef4da
# ╠═429cf33f-4422-44f0-beb8-5a1908a72273
# ╠═a963a19c-1a1f-4f44-975b-7803a6b9d9cd
# ╠═9ffe9b99-ab18-49b2-9179-3afb81b7be48
# ╟─5a212007-c0e8-4b1b-94d1-30bdb1efdb9c
# ╟─3fbc6b45-974e-430e-a4e6-960323015e74
# ╟─854731c1-7a34-4066-aa74-01629c87d75d
# ╟─13f01881-2645-429b-9856-6c3f19c0ad48
# ╟─38a45961-0ffb-43d4-aa24-36d503ed4618
# ╟─1467b184-22ac-4038-ad1b-f084d4443b27
# ╟─c87a830a-f48a-4444-81bc-3efd69a130ad
# ╟─768535e0-a514-4dff-ac8b-0d7ca126149c
# ╟─5d6222cf-99f3-4ce9-a4a2-91c17dc9c0d2
# ╟─6d37916c-7895-49d3-b8a3-c8661050ebcb
# ╟─6482d05d-06e2-43cc-ab53-ff4bbcd63e3e
# ╟─c67290fc-6291-4f3e-a660-a3c4afa3a5e3
# ╟─4e1a5050-59b0-4d24-98bb-1520c06b28c5
# ╟─42a7b186-aa04-4249-a129-bf925f181008
# ╟─e1ccfd33-3d54-4249-86f1-381a1ef90615
# ╟─b096c4f2-9dce-409d-874a-a851f577bf92
# ╟─39e50296-21ff-4407-894f-2a380dc51e21
# ╟─ca05bd4f-5656-4531-b357-331c62661174
# ╟─71eb240a-5a45-4bf3-b35c-a5820ca6da6c
# ╟─d45f106d-032a-4102-a26f-7393c2220f72
# ╟─e9216d7a-c2f3-44c0-a7d9-2c62ac35ecd9
# ╟─b4ce12e3-29ec-41ac-89d3-06d08ef2beca
