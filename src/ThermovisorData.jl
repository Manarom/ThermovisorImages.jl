module ThermovisorData
    using Images,ImageShow,ImageIO
    using Plots,CSV
    using Colors, ColorVectorSpace
    using Dates,Statistics,LinearAlgebra
    #using ImageSegmentation,IndirectArrays
    using Optim
    using LaTeXStrings
    using Distributions # to evaluate the Students coefficient
    using PerceptualColourMaps # to draw heatmaps
    using StaticArrays
    using Interpolations
    using  FileTypes
    import ImageDraw

    export RescaledImage,FilteredImage, MarkeredImage,
        full_image_flag,reduced_image_flag,draw!,
        draw,fit_centred_obj!,
        radius,diameter,area,side,center,
        CircleObj,SquareObj,RectangleObj,
        CentredObj,obj_from_vect,copyobj,
        fill_im!,fill_im_external!,
        filtered_mean,filtered_std,
        filter_image,filter_image!,
        marker_image,
        read_temperature_file,
        find_temperature_files,
        along_line_distribution,
        within_mask_line_distribution

    """
    ThermovisorData

is a package designed to process thermal images stored as matrices
Each  element of thermal image represents a temperature value. The package enables users to 
load images from files, calculate temperature distributions, and compute statistical analyses
for temperatures along specified lines. It also calculates averaged angular and radial temperature
distributions (along with standard deviations) within Regions of Interest (ROIs [`CentredObj`](@ref)) 
such as  circles, squares, and rectangles. These ROI objects can be fitted to 
distinct areas (relative to their surroundings), such as the most heated regions within
the scene.
    
    """    
    ThermovisorData

    const default_images_folder = Ref(joinpath(abspath(joinpath(@__DIR__, "..")),"thermal images"))
    const FlagMatrix = Union{Matrix{Bool},BitMatrix}
    const FlagVector = Union{Vector{Bool},BitVector}
    const int_floor_abs = Int ∘ floor ∘ abs
    const int_floor = Int ∘ floor
    const int_floor_fld = Int ∘ floor ∘ fld
    const DefColorScheme = Ref("HEAT")
    const DEFAULT_FITTING_OPTIONS = Ref(Optim.Options(x_abstol=1,iterations=30))
    include("CentredObj.jl")
    include("MarkeredImage.jl")
    """
	    RescaledImage - structure stores the image data mapped to region  [0,1]
Fields:

initial  - initial image before rescaling

sz - size of the image

min - minimum value

max  - maximum value

im - image with all values from 0 to 1 
"""
    mutable struct RescaledImage{T} 
        initial::Matrix{T}
        sz::Tuple{Int,Int}
        min::T
        max::T
        im::Matrix{T}
        RescaledImage(image::Matrix{T};negate::Bool=false) where T<:Number =begin
            sz = size(image)
            new{T}(image,sz,rescale!(copy(image))...)
        end
    end
    Base.size(image::RescaledImage) = image.sz
    Base.copy(image::RescaledImage) = RescaledImage(copy(image.initial))
    function rescale!(image::AbstractMatrix;negate::Bool=false)
        min,max = extrema(image)
        if negate 
            @. image = (image - min)/(max-min) 
        else 
            @. image = (image - min)/(max-min)
        end
        return (min,max,image)
    end
    """
        Type to store image with filtered temperature region 

:full - filtered rescaled image of the same size as the input with all pixels which are not the part of the pattern with label value 

:region_indices - cartesian indices of the pattern in the input image

:reduced - image of reduced size where all not-inpatter pixels removed  
    (the 	scaling of this image is the same as of the input `imag.initial` 
    see [`RescaledImage`](@ref) type )

:reduced_flag - bitmatrix version of (reduced)
"""
    mutable struct FilteredImage{T}
        full::RescaledImage{T}
        region_indices::Vector{CartesianIndex{2}}
        reduced::SubArray{T,2,Matrix{T},Tuple{UnitRange{Int},UnitRange{Int}},false}
        reduced_flag::SubArray{Bool,2,BitMatrix,Tuple{UnitRange{Int},UnitRange{Int}},false}
    end

    """
    full_image_flag(filtered_im::FilteredImage)

Returns the BitMatrix flag of filtered pattern in the whole image.

Can be used as index matrix in the full image e.g.:
 `filtered_image.full.initial[full_image_flag(filtered_image)]` will return 
 all elements which belong to the pattern

"""
function full_image_flag(filtered_im::FilteredImage) 
        flag = BitMatrix(undef,filtered_im.full.sz...)
        fill!(flag,false)
        @. flag[filtered_im.region_indices]=true
        return flag
    end   
    reduced_image_flag(fim::FilteredImage) =  copy(fim.reduced_flag)  
    reduced_image(fim::FilteredImage) = copy(fim.reduced)
    """
    image_discr(im1,im2)

Calculates the scalar distance between two matrices by checking the equality of their elements
"""
function image_discr(im1,im2)
        # calculates distance between two bit-images of the same size 
        N = prod(size(im1))
        return sum(1 - i[1]==i[2] for i in zip(im1,im2))/(2*N)
    end
    
    """
    draw!(image::Matrix{Float64},c::CentredObj;fill=false,thickness::Int=55,color::RGB{Float64}=RGB{Float64}(0,1,0), kwargs...)

Draws CentreObj inside the image.

image - image

c - object 

fill - if true the interior of the object will be filled 

thickness - the thickness of the object's frame

color - frame and filling color 
"""
function draw!(image::Matrix{Float64},c::CentredObj;fill=false,thickness::Int=-1,
                                        color::RGB{Float64}=RGB{Float64}(0,1,0), 
                                        color_scheme::String="",show_cross=true,kwargs...) 

        rgbim = to_rgb(image,color_scheme=color_scheme)
        #im_pic = ImageDraw.draw!(rgbim,LineTwoPoints(points_inds...), RGB{Float64}(1,0,0))             
        return   draw!(rgbim,c;
            fill=fill,thickness=thickness, 
            color=color, show_cross=show_cross,kwargs...)   
    end
    draw(image::Matrix{Float64};color_scheme::String=DefColorScheme[]) = to_rgb(image,color_scheme=color_scheme)
    draw(image::RescaledImage;color_scheme::String=DefColorScheme[]) = draw(image.initial,color_scheme=color_scheme)
    draw(image::FilteredImage;color_scheme::String=DefColorScheme[],draw_reduced::Bool=false) = draw_reduced ? draw(reduced_image(image),color_scheme=color_scheme) : draw(image.full,color_scheme=color_scheme) 
    
    function draw!(rgbim::Matrix{RGB{Float64}},
        c::CentredObj;fill=false,
        thickness::Int=55,
        color::RGB{Float64}=RGB{Float64}(0,1,0), show_cross=true,kwargs...) 

        ImageDraw.draw!(rgbim,convert_to_drawable(c,fill=fill,thickness=thickness), color; kwargs...)
        
        show_cross ? ImageDraw.draw!(rgbim, ImageDraw.Cross(ImageDraw.Point(revcentre(c)...), 50), color) : nothing

        return rgbim
    end

    """
    to_rgb(image::Matrix{Float64};color_scheme::String="")

Converts matrix to rgb martix by applyting the color scheme 
using `applycolourmap` function from `PerceptualColourMaps`  
"""
function to_rgb(image::Matrix{Float64};color_scheme::String="")
        if length(color_scheme) == 0 
            color_scheme = DefColorScheme[]
        end
        rgbimg_3D = applycolourmap(image,PerceptualColourMaps.cmap(color_scheme))

        return collect(colorview(RGB, permuteddimsview(rgbimg_3D,(3,1,2))) )
    end

    """
    draw(c::CentredObj;kwargs...)

Returns `CentredObj` image of minimal possible size
"""
function draw(c::CentredObj;kwargs...) 
        (x_left,y_left,x_right,y_right) = abs.(diagonal_points(c))
        image = fill(0.0,[y_right+y_left,x_right + x_left]...)
        image[1,1]=0.001            
        return draw!(image,c;kwargs...)
    end
    function draw_line_within_mask(image::Matrix{Float64},c::CentredObj,ang,length;thickness::Int=55,
                        color::RGB{Float64}=RGB{Float64}(0,1,0), color_scheme::String="",kwargs...) 
        rgbim = to_rgb(image,color_scheme=color_scheme)
        return draw_line_within_mask!(rgbim,c,ang,length;thickness=thickness,
                    color=color, kwargs...)
    end
    function draw_line_within_mask!(rgbim::Matrix{T},
            c::CentredObj,ang,length;
            thickness::Int=55,
            color::T=T(0,1,0), kwargs...) where T<:RGB{Float64}

            line_coords = line_within_mask(c,ang,length)
            #we should interchange the order of coordinates according to the ImageDraw demands
            ImageDraw.draw!(rgbim,ImageDraw.LineSegment(line_coords[2],line_coords[1],line_coords[4],line_coords[3]), color)
            
            return rgbim
    end

    #
    """
    fill_im_external!(img::FlagMatrix,c::CentredObj)
	
Fills image matrix `img` in a way that all pixels which are 
not within the CentreObj set to true.  See also `is_within`
"""
function fill_im_external!(img::FlagMatrix,c::CentredObj)
    for i in keys(img)
        inds = [k for k in Tuple.(i)]
        img[i] = !is_within(c,inds)
    end
    return img
end    
    """
    fit_centred_obj!(c::CentredObj,im_bin::FlagMatrix;
                                starting_point::Union{Nothing,Vector{Float64}}=nothing,
                                optimizer::Optim.ZerothOrderOptimizer = Optim.NelderMead(), 
                                options::Optim.Options=Optim.Options())

Fits [`CentredObj`](@ref) to binary image pattern (region of units) by adjusting centre coordinates and dimensions
using zeroth-order optimizers from `Optim.jl` package.

Input variables:

c - `CentredObj` (modified)   

im_bin - binarized image (BitMatrix or Matrix{Bool})

(optional)

starting_point - staring vector (uses [`fill_x0!`](@ref) function to fill starting point by default)

optimizer - zeroth-order optimizer from `Optim.jl` package

options  - optimization options from `Optim.jl` package

"""
function fit_centred_obj!(c::CentredObj,im_bin::FlagMatrix;
                                starting_point::Union{Nothing,Vector{Float64}}=nothing,
                                optimizer::Optim.ZerothOrderOptimizer = Optim.NelderMead(), 
                                options::Optim.Options=DEFAULT_FITTING_OPTIONS[]) 

        optim_fun = image_fill_discr(im_bin,c) 
        x0 = Vector{Float64}(undef,length(c))
	    if isnothing(starting_point)
		    fill_x0!(x0,im_bin,c)
        else
            x0 = copyto!(x0,starting_point)
	    end
	    optim_out = optimize(optim_fun,x0,optimizer,options)
	    return (c,Optim.minimum(optim_out),optim_out)
    end
    """
    fit_centred_obj!(c::CentredObj,image::FilteredImage;
                                        starting_point::Union{Nothing,Vector{Float64}}=nothing,fit_reduced::Bool=true,
                                        optimizer::Optim.ZerothOrderOptimizer = NelderMead(),options::Optim.Options=DEFAULT_FITTING_OPTIONS[])

Fits [`CentredObj`](@ref) (modified) to filtered image (not modified)
`fit_reduced` flag (default=true) indicates what version of the image should be fitted if true - 
reduced otherwise - full image. For other input arguments see [`fit_centred_obj!(c::CentredObj,im_bin::FlagMatrix)`](@ref)
"""
function fit_centred_obj!(c::CentredObj,image::FilteredImage;
                                        starting_point::Union{Nothing,Vector{Float64}}=nothing,fit_reduced::Bool=true,
                                        optimizer::Optim.ZerothOrderOptimizer = NelderMead(),options::Optim.Options=DEFAULT_FITTING_OPTIONS[]) 
        return fit_reduced ? fit_centred_obj!(c,reduced_image_flag(image),starting_point=starting_point,optimizer = optimizer,options=options) : fit_centred_obj!(c,full_image_flag(image),starting_point=starting_point,optimizer = optimizer,options=options)
    end
    """
    fit_all_patterns(img::RescaledImage,::Type{T}=CircleObj;
                                            level_threshold::Float64=-1.0,
                                            distance_threshold::Float64=-15.0,
                                            max_centred_objs::Int=200,
                                            sort_by_area::Bool = false,
                                            is_descend::Bool = true,
                                            optimizer::Optim.ZerothOrderOptimizer = NelderMead(),
                                            options::Optim.Options=DEFAULT_FITTING_OPTIONS[]) where T<:CentredObj

Function fits all patterns of the image `img` to the vector of [`CentredObj`](@ref) ROI objects. 
The type of ROI should be provided as a second arguments (by default it is a [`CircleObj`](@ref))

img - input image of [`RescaledImage`](@ref) type

For other input arguments see [`marker_image`](@ref) and [`fit_centred_obj!(c::CentredObj,im_bin::FlagMatrix)`](@ref)
"""
function fit_all_patterns(img::RescaledImage,::Type{T}=CircleObj;
                                            level_threshold::Float64=-1.0,
                                            distance_threshold::Float64=-15.0,
                                            max_centred_objs::Int=200,
                                            sort_by_area::Bool = false,
                                            is_descend::Bool = true,
                                            optimizer::Optim.ZerothOrderOptimizer = NelderMead(),
                                            options::Optim.Options=DEFAULT_FITTING_OPTIONS[]) where T<:CentredObj
                                            
            markers = marker_image(img,level_threshold=level_threshold,
                                    distance_threshold=distance_threshold)     

           # markers_number = count_separate_patterns(markers)       
            markers_number = length(markers) 
            is_reduced_markers_number = max_centred_objs < markers_number && max_centred_objs > 0
            if is_reduced_markers_number 
                markers_number = max_centred_objs
                if sort_by_area
                    sort_reduce!(markers,total_number=markers_number,descending=is_descend)  #sort_reduce!(m,total_number=-1,descending=descending)
                end
            elseif sort_by_area
                sort_markers!(markers,is_descend)
            end    
            if markers_number<=0 
                 return Vector{T}([])  
            else
                centered_objs_to_fit = [T() for _ in 1:markers_number] # creating empty array of centred_objects
            end
            Threads.@sync for (i,c) in enumerate(centered_objs_to_fit)
                Threads.@spawn begin 
                    (fl,i_min,i_max)  = shrinked_flag(markers,i) # returns flag and minimu and maximal indices in the initial array
                    ThermovisorData.fit_centred_obj!(c,fl,optimizer = optimizer,options = options)
                    shift!(c,i_min)
                end
            end
            return centered_objs_to_fit
    end  


    """
        `filter_image(imag::RescaledImage,markers;label=0)`

Funtion zeroes all pixels of the image, except those belonging to the specified pattern.
`image` - rescaled image (see [`RescaledImage`](@ref) type)
`markers` - the matrix of the same size as the input image, each element of this matrix has unique value-label associated with some pattern.  Function `label_components` returns the markers matrix.
(optional) - the value of the label to be selected as a pattern marker

Function returns [`FilteredImage`](@ref) object
"""
    function filter_image(imag::RescaledImage,markers::MarkeredImage;label::Int=1)
        # function extracts from markerd image
        # markers - the matrix with all labeled elements
        # imag - initial rescaled image 
        return filter_image!(copy(imag.initial),external_flag(markers,label))
    end
    """
    filter_image(imag::AbstractMatrix,c::CentredObj;external=false)

Filters image according to centered object creating new image
if external  is true than as a filtering flag the inverse of centered object image is taken
"""
   filter_image(imag::AbstractMatrix,c::CentredObj;external=false) = filter_image!(copy(imag),cent_to_flag(c,size(imag),external= !external))
   filter_image(imag,flag::FlagMatrix) =  filter_image!(copy(imag),flag)
   filter_image(imag::RescaledImage;label = 1) = filter_image(imag,marker_image(imag),label=label)
   """
    filter_image(imag::RescaledImage,c::CentredObj;external=false)

Filters image according 
"""
    filter_image(imag::RescaledImage,c::CentredObj;external=false) = filter_image(imag.initial,c;external=external)
    """
    filter_image!(imag::AbstractMatrix,flag::BitMatrix)

Returns `FilteredImage` taking all elements of imag which are not external_region_flag
"""
function filter_image!(imag::AbstractMatrix,external_region_flag::FlagMatrix)
        external_region = @view imag[external_region_flag]
        @. external_region=0.0
        @. external_region_flag = !external_region_flag
        region_area_indices = findall(external_region_flag)
        min_ind,max_ind = extrema(region_area_indices)
        square_view = @view imag[min_ind[1]:max_ind[1],min_ind[2]:max_ind[2]]
        square_view_flag =@view  external_region_flag[min_ind[1]:max_ind[1],min_ind[2]:max_ind[2]]
        return FilteredImage(RescaledImage(imag),
                region_area_indices,
                square_view, 
                square_view_flag)
    end
    filter_image!(imag::AbstractMatrix,c::CentredObj;external=false) = filter_image!(imag,cent_to_flag(c,size(imag),external=!external))
    """
    filter_image!(imag::RescaledImage{Float64},external_region_flag::FlagMatrix)::FilteredImage

In-place filtering of `RescaledImage`
"""
    function filter_image!(imag::RescaledImage{Float64},
        external_region_flag::FlagMatrix)::FilteredImage
        
        external_region = view(imag.initial,external_region_flag)
        @. external_region=0.0
        @. external_region_flag = !external_region_flag
        region_area_indices = findall(external_region_flag)
        min_ind,max_ind = extrema(region_area_indices)
        @. imag.im = imag.initial
        rescale!(imag.im)
        square_view = @view  imag.initial[min_ind[1]:max_ind[1],min_ind[2]:max_ind[2]]
        square_view_flag =@view external_region_flag[min_ind[1]:max_ind[1],min_ind[2]:max_ind[2]]
        return FilteredImage(imag,
                region_area_indices,
                square_view, 
                square_view_flag)       
    end
    filter_image!(imag::RescaledImage,c::CentredObj;external::Bool=false) = filter_image!(imag,cent_to_flag(c,imag.sz,external= !external))
    filtered_mean(fltrd::FilteredImage) = Statistics.mean(fltrd.reduced[fltrd.reduced_flag])
    filtered_std(fltrd::FilteredImage) = Statistics.std(fltrd.reduced[fltrd.reduced_flag])
    """
    marker_image(rescaled::RescaledImage,level_threshold::Float64,distance_threshold::Float64=1e-3)

Markers image patterns, input umage is `RescaledImage` image type, 
level_threshold  - should be between 0.0 and 1.0
distance_threshold  - criterium of image binarization after distance transform

returns `markers`  - matrix of Int's with the same size as the input matrix, each element 
of `markers` is the label index of individual patterns of the initial image
"""
function marker_image(rescaled::RescaledImage;
                level_threshold::Float64=-1.0,
                distance_threshold::Float64=-15.0)

    if level_threshold>1 || level_threshold <=0 # if threshold is not settled explicitly the otsu thresholding algorithm is used
        level_threshold = otsu_threshold(rescaled.im)
    end
    # thermal image is an image with several region of higher temperature
    # we want to implement the watershed algorithm to separate patterns from each other
    # first we need to negate image 
    binarized = rescaled.im .< level_threshold
    dist = distance_transform(feature_transform(binarized))
    @. dist = 1 - dist
    segments = watershed(dist, label_components(dist .< distance_threshold))

    return MarkeredImage(labels_map(segments) .* (1 .-binarized))
end


    """
    read_temperature_file(f_name::AbstractString)

Reads temeprature file `f_name` is a full file name
"""
function read_temperature_file(f_name::AbstractString)
		if isfile(f_name)
			file_full_path = f_name
		else
			file_full_path = joinpath(default_images_folder[],f_name)
            isfile(file_full_path) ? nothing : return nothing
		end
        if Is(FileType.Image, file_full_path)
            
            im_float = Float64.(Gray.(FileIO.load(file_full_path)))
        else
            im_float = CSV.Tables.matrix(CSV.File(file_full_path,header=false,types=Float64))
        end
		pic = RescaledImage(im_float);
		creation_time = mtime(file_full_path)
		return (pic,creation_time)
    end
    
    """
        `find_temperature_files(folder::AbstractString)`

Searchs the folder for thermal images files using `is_temperature_file`
Returns dictionary `Dict{String,Pair{Float64,String}}` with keys parts of files matched 
using `is_temperature_file`, values - are temperature pairs of `Float64` => `full-file-name`
When file name contains "_BB_" it supposed to be the blackbody themperature distribution       
"""
function find_temperature_files(folder::AbstractString=default_images_folder[])

        files = Dict{String,Pair{Float64,String}}()
        for file in readdir(folder)
            if !contains(file,".csv")
                continue
            end
            reg_match = is_temperature_file(file)
            if !isnothing(reg_match)
                t = parse(Float64,reg_match[1])
                t_key = contains(file,"_BB_") ? "B"*reg_match[1] : reg_match[1]
                counter = 1
                t_key_check = t_key
                while haskey(files,t_key_check)
                    t_key_check = t_key*"-"*string(counter)
                    counter+=1
                end
                t_key = t_key_check
                files[t_key] =t=>joinpath(folder,file)
            end
        end
        return files
    end
    """
        `is_temperature_file(file_name::AbstractString)`

Checks if the file with `file_name` has an appropriate name for thermovisor temperature distribution file
"""
is_temperature_file(file_name::AbstractString)=match(r"_T([1-9]|[1-9][0-9]|[1-9][0-9][0-9]|[1-9][0-9][0-9][0-9]).csv",file_name)
    
    """
    plot_along_line_distribution(along_line_length,along_line_distribution;
                                        length_scaler::Float64=1.0,
                                        is_centered::Bool=true,kwargs...)

Plots temperature distribution along the line `along_line_length` - coordinates,
`along_line_distribution` - values of temperature, `length_scaler` - length scaler 
(can be used to convert pixels to the actual length)
`is_centered` - the line length is converted to the coordinates with zero value in 
the centre of the `CentredObj`

"""
function plot_along_line_distribution(along_line_length,along_line_distribution;
                                        length_scaler::Float64=1.0,
                                        is_centered::Bool=true,kwargs...)
        #centr = center(c)
	    if !is_centered
            p_line=plot(length_scaler*along_line_length,along_line_distribution,gridlinewidth=2,framestyle = :box,kwargs...)
        else
            p_line=plot(length_scaler*(along_line_length .- along_line_length[end]/2),along_line_distribution,
            gridlinewidth=2,
            framestyle = :box,kwargs...)
        end
	    xlabel!(p_line,L"Distance \ along \ the \ line \ ,\ mm")
	    ylabel!(p_line,L"Temperature \ \degree C")
	    title!(p_line,"Temperature distribution along the line")
        return p_line
    end
    """
    plot_radial_distribution_statistics(L,mean_D::T,std_D::T,
        lower_bound::Union{T,Nothing}=nothing,upper_bound::Union{T,Nothing}=nothing;
                length_scaler::Float64=1.0,
                is_centered::Bool=true,label=nothing,
                minorgrid=true,gridlinewidth=2,title="Average temperature radial distribution",
                framestyle = :box,
                dpi=600,xlabel = "Distance  across the sample ,mm", ylabel="Temperature  °C",
                kwargs...)      where T<:AbstractVector

    Plots radial ditribution averaged value, confidence bounds and confidence bounds
    multiplied by the Student's coefficient
"""
function plot_radial_distribution_statistics(L,mean_D::T,std_D::T,
        lower_bound::Union{T,Nothing}=nothing,upper_bound::Union{T,Nothing}=nothing;
                length_scaler::Float64=1.0,
                is_centered::Bool=true,label=nothing,
                minorgrid=true,gridlinewidth=2,title="Average temperature radial distribution",framestyle = :box,
                dpi=600,xlabel = L"Distance  \ across \ the \ sample ,mm", ylabel=L"Temperature \ \degree C",
                kwargs...)      where T<:AbstractVector
        points_number = Base.length(L)
        if is_centered || length_scaler != 1.0
            L2plot = copy(L)
            if is_centered
                l_center = L[int_floor(points_number/2)] 
                @. L2plot= L-l_center
            end
            L2plot .*=length_scaler
        else
            L2plot=L
        end    
	    p=plot(L2plot,
		    mean_D,label=label,
		    minorgrid=minorgrid,
		    gridlinewidth=gridlinewidth,
		    title=title,
		    ribbon = (std_D,std_D), framestyle = framestyle,dpi=dpi,kwargs...)
	    xlabel!(p,xlabel)
	    ylabel!(p,ylabel)
        !isnothing(lower_bound) ? plot!(p,L2plot,lower_bound,linecolor=:red,label=nothing) : nothing
        !isnothing(upper_bound) ?  plot!(p,L2plot,upper_bound,linecolor=:red,label=nothing) : nothing
        return p

    end
    """
    plot_angular_distribution_statistics(angles,mean_D::T,std_D::T,
                lower_bound::Union{T,Nothing}=nothing,upper_bound::Union{T,Nothing}=nothing;
                length_scaler::Float64=1.0,
                label=nothing,
                minorgrid=true,
                gridlinewidth=2,
                title="Average temperature angular distribution",framestyle = :box,
                dpi=600,xlabel = L"Angle  ,°", ylabel=L"Temperature °C",
                kwargs...)      where T<:AbstractVector

 The same as `plot_radial_distribution_statistics` but plots averaged angular distribution
"""
function plot_angular_distribution_statistics(angles,mean_D::T,std_D::T,
                lower_bound::Union{T,Nothing}=nothing,upper_bound::Union{T,Nothing}=nothing;
                length_scaler::Float64=1.0,
                label=nothing,
                minorgrid=true,
                gridlinewidth=2,
                title="Average temperature angular distribution",framestyle = :box,
                dpi=600,xlabel = L"Angle  \ ,\degree", ylabel=L"Temperature \ \degree C",
                kwargs...)      where T<:AbstractVector

                return plot_radial_distribution_statistics(angles,mean_D,std_D,
                        lower_bound,upper_bound;
                        length_scaler=length_scaler,
                        is_centered=false,
                        label=label,
                        minorgrid=minorgrid,gridlinewidth=gridlinewidth,
                        title=title,framestyle = framestyle,
                        dpi=dpi,xlabel = xlabel, ylabel=ylabel,
                        kwargs...) 

    end
    """
    line_points_to_along_length(along_line_points::Vector{T},line_points) where T

    Converts Cartesian indices of `along_line_points` to the length along line
"""
function line_points_to_along_length(along_line_points::Vector{T},line_points) where T
        line_start = T(Tuple(line_points[1:2]))
        length_along_line = Vector{Float64}(undef,Base.length(along_line_points))
        for (i,x) in enumerate(along_line_points)
			length_along_line[i]= sqrt(sum(abs2, Tuple(x-line_start)))
	    end
        return length_along_line
    end

"""
    student_coefficient(degrees_of_freedom::Int, probability; digits::Int = 3, side::Int = 2)

    Evaluates Student's distribution coefficient
"""
function student_coefficient(degrees_of_freedom::Int, probability; digits::Int = 3, side::Int = 2)
	# dof - degrees of freedome
	# probability
	if side == 2
        probability = (1+probability)/2
    end
	return round(Distributions.quantile(Distributions.TDist(degrees_of_freedom), probability),digits=digits)
end
#
end
