### A Pluto.jl notebook ###
# v1.0.1

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

# ╔═╡ 7ea37a7f-d7df-4b82-a51e-1f0197d02213
begin 
	import Pkg
	Pkg.activate(".")
	using PlutoUI

end

# ╔═╡ 63510e0f-fb96-4e91-b5bb-44c6f54cad03
using ThermovisorImages

# ╔═╡ f6c1be87-94d2-4b08-a52d-6eb637192ee8
using Plots

# ╔═╡ b3aaa84b-215b-4822-a22c-0f8c22ce9cef
using Optim

# ╔═╡ 4460f260-f65f-446d-802c-f2197f4d6b27
md"""
### `ThermovisorImages.jl`
---------------------
**ThermovisorImages.jl** is designed to process static thermal images stored as matrices in CSV files or as image files. It treats each matrix element as a temperature value. ThermovisorImages.jl provides functions to calculate temperature distributions and perform statistical analyses of temperatures within Regions of Interest (ROIs), such as circles, squares, rectangles, or along lines. ROI objects can be fitted to image patterns (regions that stand out from the background). It is also possible to evaluate statistics across multiple ROIs, including distributions of side length, area, and perimeter.

**ThermovisorImages.jl** also provides functions to recalculate the temperature distribution of the entire image (or its part within the ROI or labeled pattern), taking into account the emissivity of the surface and the spectral range of the infrared camera.

"""

# ╔═╡ 79b31b84-afe0-4aac-90bf-97e8cbfff5e2
PlutoUI.TableOfContents()

# ╔═╡ 2c5e6e4c-92af-4991-842a-7e5bdc55a46d
md"""
###### This notebooks describes the following operations:

*  1. Loading image
*  2. Creating `CentredObj` marker (ROI object)
*  3. Finding patterns
*  4. Filtering the image
*  5. Fitting marker position and size to the pattern
*  6. Evaluating radial and angular temperature distributions
*  7. Fitting multiple ROIs to image with many patterns
*  8. Fitted ROIs spatial properties statistics
*  9. Recalculation of the temperature field of the surface and the area of interest with a new value of emissivity
"""

# ╔═╡ fc6af4b0-1127-11f0-1b66-a59d87c5b141
begin# we need the project structrue inside the notebook
	notebook_path = @__DIR__()# this notebook local path
	project_path = abspath(joinpath(notebook_path,".."))#project's path
	sources_path = joinpath(project_path,"src")# it is supposed that sources are in
	assets_folder = joinpath(project_path,"assets")
	images_folder = joinpath(project_path,"thermal images")
end;

# ╔═╡ 051044c5-760c-4b60-90fc-82a347c3b6bc
#using Revise,PlutoUI,LaTeXStrings,Images,ImageShow,Plots,BenchmarkTools,Dates,FileIO,ImageIO,Optim,CSV,Colors,ColorVectorSpace,Distributions,ColorSchemes,StaticArrays,Interpolations,FileTypes,ImageDraw,StatsBase,PlanckFunctions,RecipesBase

# ╔═╡ 4f93b7ba-3488-446d-8043-718fbdc5b808
#import Gtk

# ╔═╡ 215ed2f4-71ba-4cb5-b198-677d0d7ffb38
md" default image saving folder $(@bind image_save_folder PlutoUI.TextField(default = assets_folder))"

# ╔═╡ 870113c3-b439-4d34-90d8-fdd8a158f9dd
md"""
### 1. Loading image


All examples of temperature distribution images from real application are in 
`"...project_folder\thermal images\ "` folder. All thermal images are in csv - format

Function 
```
	ThermovisorImages.find_temperature_files(images_folder)
``` 
returns `Dict` with thermovisor data file labels matched to the full file names and actual measurements temperatures. This works simple by parsing the name of the file, to be parsed as thermal image file it's name should be like :
"any-name-T567.csv" . This file is automatically interpreted as a thermal image recorded for the temperature T = 567, here all numbers after the "T" symbol are interpreted as the temperature value.

Function 
```julia
	read_temperature_file(file_name;
							inverse_intensity::Bool=false,
							color_scheme::Symbol=:none)
``` 
can be used to read both CSV-files and images,keyword argument `color_scheme` can be used to provide the color scheme (which should be one of `ColorSheme.coloscheme` keys) to decode RGB colors to Float64 values (temperatures), by default `Gray` function is used.

"""

# ╔═╡ cd12d201-3dac-48c7-bd53-7c76944f5816
md"""
	Select this check-box to load the local file using GTK dialog box.
	$(@bind add_file PlutoUI.CheckBox(default=false))
	"""

# ╔═╡ 9fe323c0-9afc-43fd-bc21-1c45b73d50e0
md"""Set the temperature of file to be loaded $(@bind add_file_temperature PlutoUI.NumberField(0:3000,default=305))"""

# ╔═╡ dd4a9e93-0d4e-497a-8ca4-0e8f36205ffb
begin 
	files_in_dir = ThermovisorImages.find_temperature_files(images_folder)
	if add_file
		add_file_fullname = Gtk.open_dialog("select file")
		add_file_name = splitpath(add_file_fullname)[end]
		push!(files_in_dir, add_file_name=> add_file_temperature=>add_file_fullname)
	end
	
	files_in_dir_tags = keys(files_in_dir)
end

# ╔═╡ 43a1fb58-cd5e-4634-8770-0ff1809b2191
@bind temp_to_load Select([v for v in files_in_dir_tags])

# ╔═╡ 794ebd5e-e9e0-4772-98a9-43e20c7ef4da
#reading selected file
begin 
	rescaled_image = ThermovisorImages.read_temperature_file(files_in_dir[temp_to_load][2])
	rgb_image_initial = ThermovisorImages.draw(rescaled_image)
end;

# ╔═╡ 429cf33f-4422-44f0-beb8-5a1908a72273
md"""
### 2.Creating `CentredObj` marker

`CentredObj` acts as a region of interest (ROI), it has independent (from the image) coordinates and size, thus it can be used to scan the image and monitor it's values. `CentredObj`s can be of different shapes viz circle, rectangle and square. The package provides simple interface to configure custom ROI objects. 

Universal constructor allows one to create `CentredObj`s:
```julia
	obj_from_vect(::Type{T},v::AbstractVector) where T<:CentredObj
```
Here, first argument is the type of object (**`CircleObj`**, **`SquareObj`** or **`RectangleObj`**), **`v`** is the vector of parameters, first two elements of `v` are centre coordinates other are dimentions, e.g. square ROI centred at (100,120) with side length 50 can be construted by calling:

```julia
	square=obj_from_vect(SquareObj,[100,120,50])
```
To convert internal image format to the drawable RGB matrices there is a **`draw`** function, it is also possible to draw the image together with the **`CentredObj`** on it.
To convert images to RGB schemes from **[`ColorShemes.jl`](https://github.com/JuliaGraphics/ColorSchemes.jl.git)** package are used.
To draw the centred object the **[`ImageDraw.jl`](https://github.com/JuliaImages/ImageDraw.jl.git)** package is used.
```julia
	img # RescaledIMage or Matrix{Float64}
	rgb_img = ThermovisorImages.draw(img,color_scheme=:inferno)# converts to RGB
	ThermovisorImages.draw!(rgb_img,SquareObj((100,120),50)) # adds the square in image
	rgb_img = ThermovisorImages.draw(img,c::CentredObj,color_scheme=:inferno)# draws both image and ROI at once
```
"""

# ╔═╡ 7f5ec486-40d7-4e7d-9ad8-4740a1b0be22
md"""
Select marker type, 

adjust it's centre coordinates and dimentions, 

use sliders to move the ROI over the image:	

marker type = $(@bind test_mask_type Select([CircleObj=>"circle",  SquareObj =>"square",RectangleObj =>"rectangle"])) \
center point row index = $(@bind test_mask_center_x Slider(1:1000,default=200,show_value=true)) \
center point column index = $(@bind test_mask_center_y Slider(1:1000,default=200,show_value=true)) \
side a= $(@bind test_side_a Slider(1:1000,default=150,show_value=true)) 
side b= $(@bind test_side_b Slider(1:1000,default=100,show_value=true))
"""

# ╔═╡ 854731c1-7a34-4066-aa74-01629c87d75d
begin

	mm_per_pixel = Ref(0.8); # this is used to calibrate the image dimentions 
	test_coord_vect = [test_mask_center_x,test_mask_center_y,test_side_a]
	if test_mask_type==RectangleObj
		append!(test_coord_vect,test_side_b)
	end
	test_centre_obj = ThermovisorImages.obj_from_vect(test_mask_type,test_coord_vect)
	ThermovisorImages.draw!(copy(rgb_image_initial),test_centre_obj,show_cross = true,fill=false)
end

# ╔═╡ 13f01881-2645-429b-9856-6c3f19c0ad48
md"""
The following block evaluates the average and standart deviation of temperature within the ROI 

Average temperature within the marker <T>=$(ThermovisorImages.mean_within_mask(rescaled_image.initial,test_centre_obj)) 

Standat deviation of average std(T) = $(ThermovisorImages.std_within_mask(rescaled_image.initial,test_centre_obj))
"""

# ╔═╡ fd30f772-f6bc-4716-982e-e9c7fd5d5e97
plot(rescaled_image.initial[test_centre_obj],ylabel="Temperature",label=nothing)

# ╔═╡ aab55f93-1f3e-4d43-b54c-4143d6a8428d
test_centre_obj

# ╔═╡ 46e42b10-1213-48f8-a614-38c1ff86566c
md"""
`CentredObj` can also be used directly for indexing

```julia
	image[centred_obj] # returns vector of all values within the centred object
	image[centred_obj] = 10.0 # sets all pixels within the centred_obj to the same 								value
```
"""

# ╔═╡ 9f55febb-b047-4b22-8575-209d45354d51
md" Export image to file $(@bind save_image CheckBox(default = false))"

# ╔═╡ 4feea216-ee48-42a3-b4ba-454f28ff690a
begin
	if save_image 		
		FileIO.save(joinpath(image_save_folder,"initial_image.png"), rgb_image_initial)
	end
end

# ╔═╡ 5a212007-c0e8-4b1b-94d1-30bdb1efdb9c
md"""
### 3,4. Finding patterns and filtering the image
After loading from file, the image is stored as **`RescaledImage`** object, a wrapper struct, which holds both the initial image and the image with temperatures, normalized to one. To find and label image features **`marker_image`** function is used.

```julia
	marker_image(rescaled::RescaledImage;
            level_threshold::Float64=-1.0,
            distance_threshold::Float64=-15.0)
```
This function performs several operations using **`Images.jl`** package:

*  1.Binarization
*  2.Distance transform
*  3.Labeling 
*  5.Watershed

By default (if level_threshold is outside of the 0...1 range) Otsu algorithm is used to binarize the image.  It returns a **`MarkeredImage`** type object, which stores in **markers** field the matrix of labels - matrix of the same size as the initial image, but with integer values. In this matrix each pattern has individual integer value. **`MarkeredImage`** also stores the coordinates of each pattern.  
"""

# ╔═╡ 1644d6a3-8f11-4ca2-8613-b9e1d4896af0
md"""
To play with binarization adjust the following options:

level threshold $(@bind level_threshold Slider(vcat(-1.0,collect(0.0:1e-2:1)),default=-1,show_value=true))

distance threshold $(@bind distance_threshold Slider(-20.0:1.0:20,default=-15.0,show_value=true))

 The following figure shows the image of markers matrix for different binarization options. 
"""

# ╔═╡ c12addac-2880-4758-b560-42db2941a77c
begin 	
	toy_markered =ThermovisorImages.marker_image(rescaled_image,level_threshold=level_threshold,distance_threshold=distance_threshold) 
	ThermovisorImages.draw(Float64.(toy_markered.markers))
end

# ╔═╡ 3ae0c3df-7b50-4b74-b802-71931731753a
md" Number of patterns $(length(toy_markered))"

# ╔═╡ dde0003d-0fe3-41da-a582-27f88a57d2c5
md"""
There are several methods for **`filter_image`** function. 
After labeling, by default it takes the last label (maximum label value), but label value could be provided externally with a corresponding keyword **label**
If the `CentredObj` is provided as a second input argument filtering just removes all elements of the initial image which are not within the **`CentredObj`**.
Some methods of filter_image function:
```julia
	filter_image(imag::RescaledImage;label::Int=1) # filters by pattern label
	filter_image(imag::RescaledImage,c::CentredObj) # filters by roi
	filter_image(imag::RescaledImage,markers::MarkeredImage;label::Int=1) # filters by specifying the label in the markers obj
```
"""

# ╔═╡ 3fbc6b45-974e-430e-a4e6-960323015e74
md""" 

Show filtered image = $(@bind filter_init_image CheckBox(default=true))

filter by $(@bind filter_by_option Select(["CentredObj", "Pattern"]))

Show image reduced to pattern/roi size $(@bind is_show_filtered_reduced CheckBox(default=true))

"""

# ╔═╡ ca5eea20-2bb3-4407-aa09-af8de2332b84
md"##### Filtered image:"

# ╔═╡ 8b6f604d-157b-42cd-a0c6-8bd5562b47ef
begin 
	if filter_init_image
		# filtering the image
		if filter_by_option=="CentredObj"
			filtered_by_obj = ThermovisorImages.filter_image(rescaled_image,test_centre_obj)
		else
			filtered_by_obj = ThermovisorImages.filter_image(rescaled_image)
		end
		h2 = ThermovisorImages.draw(filtered_by_obj,draw_reduced=is_show_filtered_reduced) 
	end
end

# ╔═╡ 38a45961-0ffb-43d4-aa24-36d503ed4618
md"Save the current heatmap $(@bind save_distr CheckBox(default=false))"

# ╔═╡ 1467b184-22ac-4038-ad1b-f084d4443b27
save_distr ? savefig(joinpath(notebook_path,"heatmap.png")) : nothing

# ╔═╡ c87a830a-f48a-4444-81bc-3efd69a130ad
md"""

###  5,6. Evaluating radial and angular temperature distributions within the ROI
------------------------
Temperature distribution along the line, which goes through the ROI's centre within the ROI, can be obtained by calling 

```julia
	(along_line_length,distrib,line_points) = along_mask_line_distribution(imag::AbstractMatrix,c::CentredObj,direction_angle=0.0,line_length=10.0;length_per_pixel=1.0,use_wu::Bool=false)
```
It return three vectors:  along line coordinate, along line values and along line points coordinates in **`image`**, **`direction_angle`** is the rotation angle of line which goes through the center of roi and  **`line_length`** is the length of this line, additionaly, keyword argument  **`length_per_pixel`** can be provided to convert the coordinates from pixels to appropriate units.

"""

# ╔═╡ 6d37916c-7895-49d3-b8a3-c8661050ebcb
md"""
There are two algorithms available to obtain the line points within the image viz **`bresenham`** (default) and **`xiaolin_wu`**. Both were taken from the **[`ImageDraw.jl`](https://github.com/JuliaImages/ImageDraw.jl.git)** package. **`xiaolin_wu`** algorithm produces two coordinates for every point along the line, and the resulting temperature for each position is calculated as the average of temperatures for these two points. 


"""


# ╔═╡ a923a10d-064c-4403-a81c-476aef925607
md"""
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

# ╔═╡ 71eb240a-5a45-4bf3-b35c-a5820ca6da6c
md" Fit the ROI to image $(@bind is_fit_mask CheckBox(false))"

# ╔═╡ 4e1a5050-59b0-4d24-98bb-1520c06b28c5
begin 
	coord_vect = [mask_center_y,mask_center_x,side_a]
	if mask_type==RectangleObj
		append!(coord_vect,side_b)
	end
	centre_obj = ThermovisorImages.obj_from_vect(mask_type,coord_vect)
	#is_make_obj_selfy ? centre_obj_image = draw(centre_obj,thickness=25) : nothing
end;

# ╔═╡ 768535e0-a514-4dff-ac8b-0d7ca126149c
# fitting the loaded image
begin 
	filtered_by_markers = filter_image(rescaled_image,marker_image(rescaled_image))
	fitted_obj = ThermovisorImages.copyobj(centre_obj)	
	if image_type=="filtered"
		image_to_show = ThermovisorImages.reduced_image(filtered_by_markers)
	elseif image_type =="not filtered"
		image_to_show = rescaled_image.initial
	else # full filtered
		image_to_show = filtered_by_markers.full.initial
	end
	if is_fit_mask 
		fit_centred_obj!(fitted_obj,filtered_by_markers, image_type=="filtered") 
		fitted_obj
		mm_per_pixel[]=sqrt((0.25π*25^2)/ThermovisorImages.area(fitted_obj))
	else
		mm_per_pixel[]=1.0
	end
end;

# ╔═╡ 5d6222cf-99f3-4ce9-a4a2-91c17dc9c0d2
fitted_obj

# ╔═╡ e1ccfd33-3d54-4249-86f1-381a1ef90615
md"""
Upper figure shows the image, ROI and inclined line which goes through ROI's center. By adjusting ROI position, the orientation and the length of this line temperature distribution of some feature can be studied. \
 	direction angle in degrees = $(@bind direction_angle Slider(0:1:360,show_value=true,default=45)) \
	line length in $(is_fit_mask ? "mm" : "pixels" )= $(@bind line_length Slider(0.1:0.1:250,show_value=true,default=100))
"""

# ╔═╡ 42a7b186-aa04-4249-a129-bf925f181008
begin
	rgb_image = ThermovisorImages.draw(image_to_show,fitted_obj,show_cross = true)
	imag = ThermovisorImages.draw_line_within_mask!(rgb_image,fitted_obj,direction_angle,line_length/mm_per_pixel[])
	centre_obj
end

# ╔═╡ 54c27e44-d34e-49bb-9ca4-2a8218f463d5
imag

# ╔═╡ b096c4f2-9dce-409d-874a-a851f577bf92
begin 
	#ThermovisorData.diag_ang(fitted_obj)
	(along_line_length,distrib,line_points) = ThermovisorImages.along_mask_line_distribution(image_to_show,fitted_obj,direction_angle,line_length,use_wu=is_use_wu,length_per_pixel=mm_per_pixel[])
	#@show length(points)
	
	pl_distrib=plot(ThermovisorImages.along_line_distribution_plot(along_line_length,distrib,is_centered = true))
end

# ╔═╡ 59f9a7f2-9601-431c-a897-543fa25c64c4
fitted_obj

# ╔═╡ 39e50296-21ff-4407-894f-2a380dc51e21
begin 
		new_line =line_length# minimum(side(fitted_obj))-5
		ang_range = 0.0:1:180 # range of angles 
		(R,D) = ThermovisorImages.radial_distribution(image_to_show,fitted_obj,ang_range,line_length=new_line,length_per_pixel=mm_per_pixel[])

		DS = ThermovisorImages.radial_distribution_statistics(R,D)

		p_radial = plot(ThermovisorImages.radial_distribution_statistics_plot(DS,show_lower_bound=true,show_upper_bound=true,probability=0.99))

			angs = collect(ang_range)
	angular_DS = ThermovisorImages.angular_distribution_statistics(angs,R,D)
	
	p_angular = plot(ThermovisorImages.angular_distribution_statistics_plot(angular_DS))
end;

# ╔═╡ 32848d3c-866b-4a6e-be07-ff6aae73d754
md"""
It is also possible to evaluate the statistics on radial and angular distribution of image points within the ROI object.
Function 
```julia
(R,D) = radial_distribution(imag::AbstractMatrix,c::CentredObj,angles_range::AbstractRange; line_length=0.0,length_per_pixel=1.0,use_wu::Bool=false)
```
Evaluates the radial distribution, of the  **`image`** points within the object **`c`** of **`CentredObj`** subtype, **`angles_range`** is the range angles, **`line_length`** is the length of the line, **`length_per_pixel`** can be used to convert the along line coordinates from pixels to some real units. Function returns a vector **`R`** of coordinates along the line, **`D`** is the matrix, each column corresponds to the values of temperature for a particular angle in the input range.

After calculating the radial distribution; statistics on radial and angular  distribution can be evaluated by calling:

```julia
	DS_radial = radial_distribution_statistics(along_length_coordinate,
        distrib)
	DS_angular = angular_distribution_statistics(angles,along_length_coordinate,distrib)
```
These functions evaluate averaged over the angle and line length, it returns an object of **`DistributionStatistics`** type, which strores averaged values, standard deviation, and sampling volume.

Confidence boundaries according to t-distribution can be evaluated by calling **`eval_bounds`** function on `
```julia
(lower_bound,upper_bound) = eval_bounds(DS::DistributionStatistics;is_use_student::Bool=true,probability::Float64=0.95)
```
This function returns the upper and the lower values, ehich define confidence boundaries for specified probability as: 

``T = <T> \pm t \cdot S(T) ``

where S(T) is the standard deviation.

There are several predefined templates(recipes, thanks to  [`RecipesBase`](https://docs.juliaplots.org/stable/recipes) ) to plot the along line, angle- and length-averaged temperature distribution. Functions  **`along_line_distribution_plot`**, **`radial_distribution_statistics_plot`** and **`angular_distribution_statistics_plot`** return the structure which is linked to the `plot` function through the recipe macro, thus, calling `plot` will automatically setup plot properties to specific predefined values.  

```julia
using Plots
plot(along_line_distribution_plot(along_line_length,along_line_distribution))
plot(radial_distribution_statistics_plot(ds::DistributionStatistics))
plot(angular_distribution_statistics_plot(ds::DistributionStatistics))
```
These functions has several optional arguments (for details see docs), the first one plots 
"""

# ╔═╡ ea232c80-261b-4dc2-8891-2b7090f36760
p_radial

# ╔═╡ 8a558860-00d8-4f87-b900-4620881ade90
p_angular

# ╔═╡ e9216d7a-c2f3-44c0-a7d9-2c62ac35ecd9
md"Save image with marker and temperature distributions $(@bind save_average_radial_distribution CheckBox(default=false))"


# ╔═╡ b4ce12e3-29ec-41ac-89d3-06d08ef2beca
begin 
	if save_average_radial_distribution  	
		FileIO.save(joinpath(assets_folder,"filtered_image_with_marker.png"),imag)
		savefig(pl_distrib,joinpath(assets_folder,"line_distrib.png"))
		savefig(p_radial,joinpath(assets_folder,"radial_distrib.png"))
		savefig(p_angular,joinpath(assets_folder,"angular_distrib.png"))
	end
end;

# ╔═╡ cc909b53-ed4d-44a1-a410-ff25533afc2d
md"""
###  7,8. Fitting multiple ROI objects to the image with several temperature features
-----------------------
This section goes through the creation of the image with multiple randomly distributed  patterns and fitting a vector of  `CentredObj` ROIs.

First, create several  ROIs, this can be done by calling special function
```julia
	generate_random_objs(::Type{T},centers_range::NTuple{2,R},obj_number::Int,dimension_range::R) where {T<:CentredObj,R<:StepRange{Int,Int}}
```
This functions returns a Vector of `obj_number` randomly distributed in space ROIs of predefined type, `centers_range` and `dimension_range` arguments allows to provide StepRanges to set the range of ROIs centres and dimensions variation. After generating ROIs they can be placed in matrix of float by indexing. The following code places ten randome rois in the matrix. 
```julia
	im = fill(0.0,100,200)
	rois = generate_random_objs(SquareObj,(20:80,1:200),10,20:1:30)
	for r in rois
		im[r]=1.0
	end
```

 calling `fit_all_patterns` function.

```julia
	fit_all_patterns(img::RescaledImage,::Type{T}=CircleObj;
    			                    max_centred_objs::Int=200,
                			        sort_by_area::Bool = false,
                        			is_descend::Bool = true)

```
 	This function fits `CentredObj`s of specified type, to all patterns and returns the vector of ROIs. Number of ROIs can be limited to `max_centred_obj` and also it is possible to sort all patterns by area before fitting by setting `sort_by_area`, if `is_descend` is true, all patterns will be sorted in area descending order. E.g. if one needs to fit four largest patterns with rectangular ROIs:

```julia
	vect_of_rois = fit_all_patterns(img,RectangleObj,max_centred_objs=4,sort_by_area = true)
```

It is also possible to precreate the Vector of empty ROIs and fit them to the image patterns using :

```julia
	vect_of_rois = [SquareObj for  _ in 1:10] # creates 10 square ROIs
	fit_all_patterns!(vect_of_rois,marker_image(img)) # fits rois to patterns
```
"""

# ╔═╡ a76e0a08-393e-472a-8df5-0650eb6a60af
md"""
Select objects on image types = $(@bind image_patterns_type Select([CircleObj=>"circle",  SquareObj =>"square",RectangleObj =>"rectangle"])) 
"""

# ╔═╡ 6e728ea6-38be-437a-96b4-9fa084f8fec5
md"""
Select objects ROI type = $(@bind multifit_roi_type Select([CircleObj=>"circle",  SquareObj =>"square",RectangleObj =>"rectangle"])) 
"""

# ╔═╡ 0badf26a-38fa-45be-9704-d4e80b12a9cb
md"""
	Select this checkbox to fit all objects $(@bind is_fit_multiple CheckBox(default=true))
	"""

# ╔═╡ f8154558-d0cb-4b27-8c0d-b5cac07a099c
md"Patterns size range $(@bind pattern_sizes_range confirm(RangeSlider(5:1:120,default=10:1:50)))"

# ╔═╡ 39e31290-e7b5-47ce-ac46-aebc33ddfa54
md"""
	Number of generated patterns $(@bind patterns_number confirm(PlutoUI.Slider(1:1:200,default=30,show_value=true)))
	"""

# ╔═╡ 7e484bde-1f2f-4d32-87c6-64ff884dc272
md""" Number of patterns to fit $(@bind max_obj_number  Select(1:patterns_number,default=patterns_number))"""

# ╔═╡ 304ac75d-bdc4-41de-bd2f-d1843ecd22f9
md"Sort by area? $(@bind is_sort_by_area CheckBox(false))"

# ╔═╡ 10954f10-9414-4839-872f-c2516d5d8e4e
md"""
	Click this button to generate the image pattern $(@bind fit_multiple Button("Regenerate patterns"))
	"""

# ╔═╡ d5b6f453-5e92-41e6-a45f-cb75660bc198
begin # generating pattern 
	fit_multiple
	#patterns_number = 100
	image_size = (500,1000)
	img = fill(0.0,image_size...)# filling initial scene
	generated_rois = ThermovisorImages.generate_random_objs(image_patterns_type,(1:1:image_size[1],1:1:image_size[2]),patterns_number,pattern_sizes_range)
	for r in generated_rois
		img[r] = rand(5:10)
	end
	rs = RescaledImage(img)
	markers = ThermovisorImages.marker_image(rs , level_threshold = 0.2)
	separate_patterns_number = length(markers) 
	# now we need to check for the separate patterns number 
	rgb_markers = ThermovisorImages.draw(rs)	
	# converting to rgb image
end;

# ╔═╡ adc6d2f9-e247-4447-8d97-67c9c1f54135
ThermovisorImages.marker_image(rs, level_threshold = 0.2)

# ╔═╡ f27fb11f-66f7-45b7-9cb3-75dc1f9b05d8
draw(rs)

# ╔═╡ 9c13d94e-ca2f-41d6-922a-428bb7a476c8
generated_patterns_stat = ThermovisorImages.CentredObjCollectionStat(generated_rois)

# ╔═╡ 96ad6e27-52dd-41aa-b115-f852049a485a
md"""Number of separate patterns:    $(separate_patterns_number)"""

# ╔═╡ 9d79ba97-5aa8-4d60-8cf4-523a28b2e5ae
md"Generated patterns"

# ╔═╡ 6c78d805-4b14-4b7f-ad96-439d2a56605e
rgb_markers

# ╔═╡ 68f1ba03-a84a-40fb-9ce4-bf0ac9ae55d0
md" Distance threshold $(@bind multifit_distance_threshold Slider(-15.0:1e-1:10.0,default=-15.0,show_value=true))"

# ╔═╡ 6adfae4d-5137-4692-b9f3-3793c4c76202
begin # fitting ROI's to image with several 
	fit_multiple
	if is_fit_multiple

		fitted_rois = ThermovisorImages.fit_all_patterns(rs, multifit_roi_type, sort_by_area=is_sort_by_area, max_centred_objs=max_obj_number, distance_threshold=multifit_distance_threshold  , level_threshold = 0.1)

	else
		println("There is no fitted ROIs")
	end
end;

# ╔═╡ fcb71c81-8ee5-4cf7-b293-ab97261d7213
md"Fitted $(length(fitted_rois)) $(is_sort_by_area ? :largest : :random) patterns"

# ╔═╡ 550a5eaf-d608-47ae-b444-10aab596ef3d
begin
if length(fitted_rois)>0
			rgb_image_multi_roi = ThermovisorImages.draw(img,fitted_rois[1],show_cross = true,fill=true)
			for i in 2:length(fitted_rois)
			 	ThermovisorImages.draw!(rgb_image_multi_roi,fitted_rois[i],fill=true,show_cross = true)
			end
			 rgb_image_multi_roi
		end
end

# ╔═╡ 0cdb27f4-9796-479b-b43f-b349eaabc049
md"""
After fitting the regions of interest, it is possible to calculate basic statistics on the geometric properties of the objects by creating the objects of `CentredObjCollectionStat` type. This object contains data on maximal, minimal,averages, standard deviations and disribution hystorgams of areas, perimeters and sides of ROIs.
```julia
	stat_on_rois = CentredObjCollectionStat(fitted_rois)
	stat_on_rois.maxs.side # returns the Tuple of (maximum side values,object index)
	stat_on_rois.maxs.area # the same, but for area
	h = stat_on_rois.hists.area # return histogram which can be plotted using plot(h)
```
and plot the histogram
```julia
	plot(stat_on_rois.hists.area) # plots side - distribution histogram 
```


"""

# ╔═╡ f1512fc5-4a7a-4274-9d42-3057d9aec04f
StatOnRois = ThermovisorImages.CentredObjCollectionStat(fitted_rois,nbins=Int(floor(patterns_number/3)))

# ╔═╡ d91b478b-57fa-4f26-a05c-649097202102
md"""
	Maximal area of  $(StatOnRois.maxs.area[1])  has the ROI # $(StatOnRois.maxs.area[2])

	The smallest ROIs # $(StatOnRois.mins.area[2]) has the area of  $(StatOnRois.mins.area[1])

	"""

# ╔═╡ 2925fafa-4722-4335-ba49-77c6a8fb110b
md"""
Show histogram for $(@bind selected_hist Select([:side,:area,:perimeter]))

"""

# ╔═╡ e7e6884a-1145-4a01-a429-6c4a84e7ea33
multi_hist = plot(getfield(StatOnRois.hists,selected_hist),title="Distribution of ROIs over  "*string(selected_hist),label=nothing,alpha=0.5,dpi=600)

# ╔═╡ 68b33b39-5ef5-4560-b4b2-1fe2f43a3628
md" Save multiple patterns fit ? $(@bind is_save_multipattern_fit CheckBox(default=false))"

# ╔═╡ 8a132aba-aa8a-428a-84a2-0ab6e5e2b891
begin 
	if is_save_multipattern_fit
		FileIO.save(joinpath(assets_folder,"multiple_patterns_initial.png"),rgb_markers)
		if is_fit_multiple
			FileIO.save(joinpath(assets_folder,"multiple_patterns_fitted.png"),rgb_image_multi_roi)
			savefig(multi_hist,joinpath(assets_folder,"multiple_patterns_hist.png"))
		end
	end
end;

# ╔═╡ f357bbd8-4d81-4f9e-870e-cf57124c5042
md"""
The following blocks show the same usage as the previous one, with an axception that it analyzes standard image files.  
"""

# ╔═╡ 821a7c95-f4da-410d-b780-111abb6d0db5
md"""
	Show/hide fitted rois $(@bind is_draw_rois Select( ["show" ;"hide" ],default = :hide))
	"""

# ╔═╡ febd591e-bb9f-4b21-93c8-aafd4c81ce12
begin
	coins_file = joinpath(images_folder,"coins.jpg")
	if !isfile(coins_file)
		download("http://docs.opencv.org/3.1.0/water_coins.jpg",coins_file)
	end	
	rescaled_coin = ThermovisorImages.read_temperature_file(coins_file,inverse_intensity=true)
	markers_coins = ThermovisorImages.marker_image(rescaled_coin)
	im_coin_float = rescaled_coin.im
	fitted_rois_coins = ThermovisorImages.fit_all_patterns(rescaled_coin,multifit_roi_type)
end;

# ╔═╡ 0044c49b-1c72-4f78-97ee-87932c97d2a9
begin
		if is_draw_rois=="show"
		rgb_image_coins = ThermovisorImages.draw(im_coin_float,fitted_rois_coins[1],show_cross = true,fill=true)
		for i in 2:length(fitted_rois_coins)
				 ThermovisorImages.draw!(rgb_image_coins,fitted_rois_coins[i],thickness = 1,fill=true,show_cross = true)
		end
	else
		rgb_image_coins = ThermovisorImages.to_rgb(im_coin_float)
		
	end
	rgb_image_coins
end

# ╔═╡ b640fcd0-3e49-471d-b281-87137a781eba
begin 
	test_markers = ThermovisorImages.marker_image(rescaled_coin)
	before_sorting = ThermovisorImages.draw(Float64.(test_markers.markers))
	ThermovisorImages.sort_markers!(test_markers)
	after_sorting = ThermovisorImages.draw(Float64.(test_markers.markers))
end;

# ╔═╡ 8b2fdcf1-cd0e-4234-a24a-afa597552f9e
md""" After labeling patterns (without fitting the ROIs), labels in `MarkeredImage` can be sorted by their area. Each pattern is assigned an integer label according to the internal logic of the labeling algorithm. In some cases, it may be useful to sort the labels by the size of the patterns, for example, before ROI fitting, when the number of ROIs is limited to a value smaller than the total number of patterns.  

```julia
markered_image = marker_image(rescaled_image)# markers are not ordered   
sort_markers!(markered_image) # now,the pattern with index `1` has the largest area 
rois = [SquareObj for _ in 1:5] # five empty squares
fit_all_patterns!(rois,markered_image) # now five largeest pattern are fitted
```
"""

# ╔═╡ 054b15d8-a4e6-42d4-b097-938d05cbb198
md"Markered image before (left) and after (right) sorting features by their area (the brighter the color, the higher the label of the pattern)" 

# ╔═╡ 7a00ce43-94e2-4f68-b651-b57bf7d6ab05
hcat(before_sorting,after_sorting)

# ╔═╡ 20cf3079-1115-451b-870d-2457a5cfd333
md""" 
### 9. Recalculation of the temperature field of the surface and the area of interest with a new value of emissivity. 
----------------------------------
An infrared camera measures the integrated radiance from a surface within its operational spectral range and converts this measurement into temperature using its calibration matrix. This measured value depends on both the surface’s emissivity and its temperature. When the emissivity is adjusted—such as through correction or updating—the temperature must be recalculated to ensure that the modeled radiance aligns with the observed radiance. This recalculation is essential because emissivity varies depending on the material, surface condition, and wavelength; using an incorrect emissivity value results in errors in temperature measurement.  

The main idea is to recalculate the temperature in such a way that the measured radiance within the infrared camera’s spectral range remains the same for two different surface emissivities (initial and modified).

To obtain the new temperature ``T_{new}``, the package solves the following nonlinear equation using univariate optimization method:

``\epsilon_{initial} \int_{\lambda_{left}}^{\lambda_{right}}I_{bb}(\lambda,T_{initial})d\lambda =\epsilon_{new} \int_{\lambda_{left}}^{\lambda_{right}}I_{bb}(\lambda,T_{new})d\lambda``

where ``\lambda_{left}`` and ``\lambda_{right}`` define the working range of thermographic camera, ``T_{initial}`` and ``\epsilon_{initial}`` are the temperature and emissivity at each point of image,  respectively, ``\epsilon_{new}`` is a new emissivity.


**ThermovisorImages.jl** package provides `recalculate_with_new_emissivity` function to accomplish this task. This function has three methods:

```julia
recalculate_with_new_emissivity!(image::AbstractArray,new_emissivity::Float64,
                                        image_emissivity::Float64;
                                        λ_left::Union{Float64,Nothing}=14.0,
                                        λ_right::Union{Float64,Nothing}=14.5,
                                        is_in_kelvins::Bool=false,
                                        rel_tol::Float64=1e-3) 				 #(1)
recalculate_with_new_emissivity!(image::AbstractArray,c::CentredObj,
									new_emissivity::Float64, image_emissivity::Float64;kwargs...)     #(2)
recalculate_with_new_emissivity!(image::AbstractArray,marker::MarkeredImage,
        								label::Int,new_emissivity::Float64,
										image_emissivity::Float64;kwargs...) #(3)
```
First method does not depend on package's internals and recalculates the whole image (or image view provided externally). `new_emissivity` and `image_emissivity` are the initial emissivity of the image and new emissivity for which the temperature should be adjusted. If both λ_left and λ_right are equal, function uses single wavelegth equation:

``\epsilon_{initial} I_{bb}(\lambda_{fixed},T_{initial}) =\epsilon_{new}I_{bb}(\lambda,T_{new})``

Singlewavelength version is much faster than the one with band integration, but less accurate.

Methods (2) and (3) both allow to recalculate the image within the specified region, (2)  - recalculates the temperature within the ROI and (3)  - within the pattern, labeled with the `label` value.

"""

# ╔═╡ 0e05f2b9-37d2-4626-b50c-4c8d48022904
md"""
Select thermal image by tag $(@bind temp_to_recalc Select([v for v in files_in_dir_tags]))
"""

# ╔═╡ dc5be80b-9a5e-42f2-b75f-b338292851ee
#reading selected file
begin 
	im_to_recalc = ThermovisorImages.read_temperature_file(files_in_dir[temp_to_recalc][2])
	roi_obj = ThermovisorImages.fit_all_patterns(im_to_recalc)
end;

# ╔═╡ da18cd4d-73b3-491f-b9f0-d374b92ed8d2
@bind args_in confirm(PlutoUI.combine() do Child
md"""
	Infrared camera settings:
	
	``\lambda_{left}`` = $(Child(Slider(0.5:1e-1:30,default=14.0,show_value=true)))
	
	``\lambda_{right}`` = $(Child(Slider(0.5:1e-1:30,default=14.5,show_value=true)))

	``new \ emissivity`` = $(Child(Slider(0.01:1e-2:1.0,default=0.70,show_value=true)))	
	
	``\epsilon_{image}`` = $(Child(Slider(0.01:1e-2:1.0,default=1.0,show_value=true)))

	Use single wavelength ? (calculating in band can be time consuming) = $(Child(CheckBox(true)))
	
	"""

end)

# ╔═╡ 8b017646-7a75-4973-a204-d74a42ffc97f
args_t = NamedTuple{(:λ_left,:λ_right,:new_emissivity,:image_emissivity)}( !args_in[5] ? args_in[1:4] :  (args_in[1],args_in[1],args_in[3:4]...) )

# ╔═╡ 4d377da2-e1b2-4aa3-a2c5-6d176ae8905f
md"Left  - initial temperature distribution, right - recalculated "

# ╔═╡ 9ba9540f-a2e2-40c4-99d5-940ec2e2839b
begin
	im_copy = copy(im_to_recalc.initial)
	ThermovisorImages.recalculate_with_new_emissivity!(im_copy,roi_obj[],args_t.new_emissivity, args_t.image_emissivity;λ_left = args_t.λ_left, λ_right = args_t.λ_right)
	heatmap(hcat(im_to_recalc.initial,im_copy),aspect_ratio=:equal,title = "Left - initial, right - recalculated")
end

# ╔═╡ 577897c4-c042-495e-a10e-9ac07ab2bf2b
md"Difference between initial temperature distribution and recaculated for the new value of emissivity"

# ╔═╡ d7abe315-ad5e-485e-8873-fe6f3cd241b5
heatmap(im_copy .-im_to_recalc.initial,title="Temperature difference")

# ╔═╡ Cell order:
# ╟─4460f260-f65f-446d-802c-f2197f4d6b27
# ╠═7ea37a7f-d7df-4b82-a51e-1f0197d02213
# ╟─79b31b84-afe0-4aac-90bf-97e8cbfff5e2
# ╟─2c5e6e4c-92af-4991-842a-7e5bdc55a46d
# ╟─fc6af4b0-1127-11f0-1b66-a59d87c5b141
# ╠═051044c5-760c-4b60-90fc-82a347c3b6bc
# ╠═4f93b7ba-3488-446d-8043-718fbdc5b808
# ╟─215ed2f4-71ba-4cb5-b198-677d0d7ffb38
# ╠═63510e0f-fb96-4e91-b5bb-44c6f54cad03
# ╠═f6c1be87-94d2-4b08-a52d-6eb637192ee8
# ╠═b3aaa84b-215b-4822-a22c-0f8c22ce9cef
# ╟─870113c3-b439-4d34-90d8-fdd8a158f9dd
# ╟─dd4a9e93-0d4e-497a-8ca4-0e8f36205ffb
# ╟─43a1fb58-cd5e-4634-8770-0ff1809b2191
# ╟─cd12d201-3dac-48c7-bd53-7c76944f5816
# ╟─9fe323c0-9afc-43fd-bc21-1c45b73d50e0
# ╠═794ebd5e-e9e0-4772-98a9-43e20c7ef4da
# ╟─429cf33f-4422-44f0-beb8-5a1908a72273
# ╟─7f5ec486-40d7-4e7d-9ad8-4740a1b0be22
# ╟─13f01881-2645-429b-9856-6c3f19c0ad48
# ╟─fd30f772-f6bc-4716-982e-e9c7fd5d5e97
# ╟─854731c1-7a34-4066-aa74-01629c87d75d
# ╟─aab55f93-1f3e-4d43-b54c-4143d6a8428d
# ╟─46e42b10-1213-48f8-a614-38c1ff86566c
# ╟─9f55febb-b047-4b22-8575-209d45354d51
# ╟─4feea216-ee48-42a3-b4ba-454f28ff690a
# ╟─5a212007-c0e8-4b1b-94d1-30bdb1efdb9c
# ╟─1644d6a3-8f11-4ca2-8613-b9e1d4896af0
# ╟─3ae0c3df-7b50-4b74-b802-71931731753a
# ╟─c12addac-2880-4758-b560-42db2941a77c
# ╟─dde0003d-0fe3-41da-a582-27f88a57d2c5
# ╟─3fbc6b45-974e-430e-a4e6-960323015e74
# ╟─ca5eea20-2bb3-4407-aa09-af8de2332b84
# ╟─8b6f604d-157b-42cd-a0c6-8bd5562b47ef
# ╟─38a45961-0ffb-43d4-aa24-36d503ed4618
# ╟─1467b184-22ac-4038-ad1b-f084d4443b27
# ╟─c87a830a-f48a-4444-81bc-3efd69a130ad
# ╟─768535e0-a514-4dff-ac8b-0d7ca126149c
# ╟─5d6222cf-99f3-4ce9-a4a2-91c17dc9c0d2
# ╟─6d37916c-7895-49d3-b8a3-c8661050ebcb
# ╟─a923a10d-064c-4403-a81c-476aef925607
# ╟─6482d05d-06e2-43cc-ab53-ff4bbcd63e3e
# ╟─c67290fc-6291-4f3e-a660-a3c4afa3a5e3
# ╟─71eb240a-5a45-4bf3-b35c-a5820ca6da6c
# ╟─4e1a5050-59b0-4d24-98bb-1520c06b28c5
# ╟─42a7b186-aa04-4249-a129-bf925f181008
# ╟─54c27e44-d34e-49bb-9ca4-2a8218f463d5
# ╟─e1ccfd33-3d54-4249-86f1-381a1ef90615
# ╟─b096c4f2-9dce-409d-874a-a851f577bf92
# ╟─59f9a7f2-9601-431c-a897-543fa25c64c4
# ╟─39e50296-21ff-4407-894f-2a380dc51e21
# ╟─32848d3c-866b-4a6e-be07-ff6aae73d754
# ╟─ea232c80-261b-4dc2-8891-2b7090f36760
# ╟─8a558860-00d8-4f87-b900-4620881ade90
# ╟─e9216d7a-c2f3-44c0-a7d9-2c62ac35ecd9
# ╟─b4ce12e3-29ec-41ac-89d3-06d08ef2beca
# ╟─cc909b53-ed4d-44a1-a410-ff25533afc2d
# ╠═d5b6f453-5e92-41e6-a45f-cb75660bc198
# ╠═adc6d2f9-e247-4447-8d97-67c9c1f54135
# ╠═f27fb11f-66f7-45b7-9cb3-75dc1f9b05d8
# ╟─a76e0a08-393e-472a-8df5-0650eb6a60af
# ╟─6e728ea6-38be-437a-96b4-9fa084f8fec5
# ╟─0badf26a-38fa-45be-9704-d4e80b12a9cb
# ╟─f8154558-d0cb-4b27-8c0d-b5cac07a099c
# ╟─39e31290-e7b5-47ce-ac46-aebc33ddfa54
# ╠═9c13d94e-ca2f-41d6-922a-428bb7a476c8
# ╠═96ad6e27-52dd-41aa-b115-f852049a485a
# ╟─7e484bde-1f2f-4d32-87c6-64ff884dc272
# ╟─304ac75d-bdc4-41de-bd2f-d1843ecd22f9
# ╟─10954f10-9414-4839-872f-c2516d5d8e4e
# ╟─9d79ba97-5aa8-4d60-8cf4-523a28b2e5ae
# ╠═6c78d805-4b14-4b7f-ad96-439d2a56605e
# ╟─fcb71c81-8ee5-4cf7-b293-ab97261d7213
# ╟─68f1ba03-a84a-40fb-9ce4-bf0ac9ae55d0
# ╟─550a5eaf-d608-47ae-b444-10aab596ef3d
# ╟─6adfae4d-5137-4692-b9f3-3793c4c76202
# ╟─0cdb27f4-9796-479b-b43f-b349eaabc049
# ╟─d91b478b-57fa-4f26-a05c-649097202102
# ╟─f1512fc5-4a7a-4274-9d42-3057d9aec04f
# ╟─2925fafa-4722-4335-ba49-77c6a8fb110b
# ╟─e7e6884a-1145-4a01-a429-6c4a84e7ea33
# ╟─68b33b39-5ef5-4560-b4b2-1fe2f43a3628
# ╟─8a132aba-aa8a-428a-84a2-0ab6e5e2b891
# ╟─f357bbd8-4d81-4f9e-870e-cf57124c5042
# ╟─821a7c95-f4da-410d-b780-111abb6d0db5
# ╠═febd591e-bb9f-4b21-93c8-aafd4c81ce12
# ╠═0044c49b-1c72-4f78-97ee-87932c97d2a9
# ╟─b640fcd0-3e49-471d-b281-87137a781eba
# ╟─8b2fdcf1-cd0e-4234-a24a-afa597552f9e
# ╟─054b15d8-a4e6-42d4-b097-938d05cbb198
# ╟─7a00ce43-94e2-4f68-b651-b57bf7d6ab05
# ╟─20cf3079-1115-451b-870d-2457a5cfd333
# ╟─0e05f2b9-37d2-4626-b50c-4c8d48022904
# ╟─dc5be80b-9a5e-42f2-b75f-b338292851ee
# ╟─da18cd4d-73b3-491f-b9f0-d374b92ed8d2
# ╟─8b017646-7a75-4973-a204-d74a42ffc97f
# ╟─4d377da2-e1b2-4aa3-a2c5-6d176ae8905f
# ╟─9ba9540f-a2e2-40c4-99d5-940ec2e2839b
# ╟─577897c4-c042-495e-a10e-9ac07ab2bf2b
# ╟─d7abe315-ad5e-485e-8873-fe6f3cd241b5
