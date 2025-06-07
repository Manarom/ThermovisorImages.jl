 
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
            markers_number>0 ? nothing : return Vector{T}([])
            is_reduced_markers_number = max_centred_objs < markers_number && max_centred_objs > 0
            if is_reduced_markers_number 
                markers_number = max_centred_objs
                if sort_by_area
                    sort_reduce!(markers,total_number=markers_number,descending=is_descend)  #sort_reduce!(m,total_number=-1,descending=descending)
                end
            elseif sort_by_area
                sort_markers!(markers,is_descend)
            end    
            centered_objs_to_fit = [T() for _ in 1:markers_number] # creating empty array of centred_objects

            Threads.@sync for (i,c) in enumerate(centered_objs_to_fit)
                Threads.@spawn begin 
                    (fl,i_min,i_max)  = shrinked_flag(markers,i) # returns flag and minimu and maximal indices in the initial array
                    fit_centred_obj!(c,fl,optimizer = optimizer,options = options)
                    shift!(c,i_min - CartesianIndex(1,1))
                end
            end
            return centered_objs_to_fit
    end  
    """
    image_discr(im1,im2)

Calculates the scalar distance between two matrices by checking the equality of their elements
"""
function image_discr(im1,im2)
        # calculates distance between two bit-images of the same size 
        N = prod(size(im1))
        return sum(1 - i[1]==i[2] for i in zip(im1,im2))/(2*N)
    end


    