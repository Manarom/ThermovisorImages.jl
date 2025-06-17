 
 export fit_centred_obj!,fit_all_patterns!,fit_all_patterns
 """
    fit_centred_obj!(c::CentredObj,im_bin::FlagMatrix;
                                starting_point::Union{Nothing,Vector{Float64}}=nothing,
                                optimizer::Optim.ZerothOrderOptimizer = Optim.NelderMead(), 
                                options::Optim.Options=DEFAULT_FITTING_OPTIONS[],refit::Bool = true)

Fits [`CentredObj`](@ref) to binary image pattern (region of units) by adjusting centre coordinates and dimensions
using zeroth-order optimizers from `Optim.jl` package.

Input variables:

c - `CentredObj` (modified)   

im_bin - binarized image (BitMatrix or Matrix{Bool})

(optional)

starting_point - staring vector (uses [`fill_x0!`](@ref) function to fill starting point by default)

optimizer - zeroth-order optimizer from `Optim.jl` package

options  - optimization options from `Optim.jl` package

refit - if true the starting point of the optimization recalculated otherwise it is taken from the current state vector of ROI

"""
function fit_centred_obj!(c::CentredObj,im_bin::FlagMatrix;
                                starting_point::Union{Nothing,Vector{Float64}}=nothing,
                                optimizer::Optim.ZerothOrderOptimizer = Optim.NelderMead(), 
                                options::Optim.Options=DEFAULT_FITTING_OPTIONS[],
                                refit::Bool = true) 

        optim_fun = image_fill_discr(im_bin,c) 
        x0 = MVector{length(c),Float64}(undef)
	    if isnothing(starting_point)
            if !refit
		        fill_x0!(x0,im_bin,c)
            else
                fill_vect!(x0, c)
            end
        else
            x0 = copyto!(x0,starting_point)
	    end
	    optim_out = optimize(optim_fun,x0,optimizer,options)
	    return (c,Optim.minimum(optim_out),optim_out)
    end
    """
    fit_centred_obj!(c::CentredObj,image::FilteredImage;kwargs...)

Fits [`CentredObj`](@ref) (modified) to filtered image (not modified)
`fit_reduced` flag (default=true) indicates what version of the image should be fitted if true - 
reduced otherwise - full image. For other input arguments see [`fit_centred_obj!`](@ref)
"""
function fit_centred_obj!(c::CentredObj,image::FilteredImage,fit_reduced::Bool=true;kwargs...) 
        return fit_reduced ? fit_centred_obj!(c,reduced_image_flag(image);kwargs...) : fit_centred_obj!(c,full_image_flag(image);kwargs...)
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

For other input arguments see [`marker_image`](@ref) and [`fit_centred_obj!`](@ref)
"""
function fit_all_patterns!(markers::MarkeredImage,::Type{T}=CircleObj;
                                            max_centred_objs::Int=200,
                                            sort_by_area::Bool = false,
                                            is_descend::Bool = true,
                                            optimizer::Optim.ZerothOrderOptimizer = NelderMead(),
                                            options::Optim.Options=DEFAULT_FITTING_OPTIONS[],
                                            refit::Bool=false) where T<:CentredObj
                                            
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
            c_vect = [T() for _ in 1:markers_number] # creating empty array of centred_objects
            return fit_all_patterns!(c_vect, markers, optimizer,options,refit)
    end  
    """
    fit_all_patterns!(c_vect::Vector{T},
                            markers::MarkeredImage,
                            optimizer::Optim.ZerothOrderOptimizer = NelderMead(),
                            options::Optim.Options=DEFAULT_FITTING_OPTIONS[]) where T<:CentredObj

Fits markered image pattern and fills precreated avector of `Centerdobj`s
See [`fit_all_patterns!`](@ref)
"""
function fit_all_patterns!(c_vect::Vector{T},
                            markers::MarkeredImage,
                            optimizer::Optim.ZerothOrderOptimizer = NelderMead(),
                            options::Optim.Options=DEFAULT_FITTING_OPTIONS[],
                            refit::Bool = false) where T<:CentredObj

            num_patterns = length(markers)
            c_length = length(c_vect) 
            num_patterns =  num_patterns < c_length ? num_patterns : c_length
            v = @view c_vect[1:num_patterns]
            Threads.@sync for (i,c) in enumerate(v)
                Threads.@spawn begin 
                    (fl,i_min,_)  = shrinked_flag(markers,i) # returns flag and minimu and maximal indices in the initial array
                    fit_centred_obj!(c, fl, optimizer = optimizer, options = options, refit = refit)
                    shift!(c,i_min - CartesianIndex(1,1))
                end
            end
            return c_vect
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

Markers `RescaledImage` and fits all patterns
See [`fit_all_patterns!`](@ref)
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

        return length(markers)>0 ? fit_all_patterns!(markers,T;
                                max_centred_objs=max_centred_objs,
                                sort_by_area=sort_by_area,
                                is_descend = is_descend,
                                optimizer=optimizer,
                                options=options) : Vector{T}([])
    end
    """
    image_discr(im1,im2)

Calculates the scalar distance between two matrices by checking the equality of their elements
"""
function image_discr(im1,im2)
        # calculates distance between two bit-images of the same size 
        N = prod(size(im1))
        d = 0
        for i in 1:length(im1)
            @inline d += im1[i]!=im2[i]
        end
        return d/(2*N)#sum(1 - i[1]==i[2] for i in zip(im1,im2))/(2*N)
    end


    """
    image_fill_discr(image::FlagMatrix,c::CentredObj)

Function returns the function to evaluate the discrepancy  between 
`CentredObj` and the matrix, this function is used during the fitting procedure 
"""    
function image_fill_discr(image::FlagMatrix,c::CentredObj)
     discr_obj =ImCentDiscr(image,c)   
     return discr_obj
                #image_discr(image, fill_im!(im_copy,fill_from_vect!(c,x)))

end
struct ImCentDiscr{T<:CentredObj,F<:FlagMatrix}
    image_target::F
    image_fillable::F
    target_area::Float64
    c::T
    ImCentDiscr(imag::F,c::T) where {T<:CentredObj,F<:FlagMatrix}=  begin
        new{T,F}(imag,copy(imag),sum(imag),c)
    end
end
(icd::ImCentDiscr{SquareObj,F})(x::AbstractVector) where F<:FlagMatrix = begin
    fill_from_vect!(icd.c,x)
    s = area(icd.c)
    fill_im!(icd.image_fillable,icd.c)
    #s_im = sum(icd.image_fillable)
    return +(abs(s-icd.target_area)/icd.target_area , image_discr(icd.image_target,icd.image_fillable))
end
(icd::ImCentDiscr)(x::AbstractVector) = begin
    fill_from_vect!(icd.c,x)
    s = area(icd.c)
    fill_im!(icd.image_fillable,icd.c)
    s_im = sum(icd.image_fillable)
    return +(abs(s-s_im)/s_im , image_discr(icd.image_target,icd.image_fillable))
end
(icd::ImCentDiscr{CircleObj,F})(x::AbstractVector) where F<:FlagMatrix = begin
    fill_im!(icd.image_fillable,fill_from_vect!(icd.c,x))
    return image_discr(icd.image_target,icd.image_fillable)
end
    