#brings custom types of images 
# RescaledImage - normalized image version
# FilteredImage - 
export RescaledImage,FilteredImage,filter_image,
filter_image!,filtered_mean,filtered_std
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
        inverted::Bool
        RescaledImage(image::Matrix{T};inverse_intensity::Bool=false) where T<:Number =begin
            sz = size(image)
            new{T}(image,sz,rescale!(copy(image),inverse_intensity=inverse_intensity)...)
        end
    end
    Base.size(image::RescaledImage) = image.sz
    Base.copy(image::RescaledImage) = RescaledImage(copy(image.initial))
    function rescale!(image::AbstractMatrix;inverse_intensity::Bool=false)
        min,max = extrema(image)
        if inverse_intensity 
            @. image =(max - image)/(max-min)
        else 
            @. image = (image - min)/(max-min) 
        end
        return (min,max,image,inverse_intensity)
    end
    """
        Type to store image with filtered temperature region 

Fields:

`full` - filtered rescaled image of the same size as the input with all pixels which are not the part of the pattern with label value 

`region_indices` - cartesian indices of the pattern in the input image

`reduced` - image of reduced size where all not-inpatter pixels removed  
    (the 	scaling of this image is the same as of the input `imag.initial` 
    see [`RescaledImage`](@ref) type )

`reduced_flag` - bitmatrix or Matrix{Bool} version of (reduced image)

"""
    mutable struct FilteredImage{T}
        full::RescaledImage{T}
        region_indices::Vector{CartesianIndex{2}}
        reduced::SubArray{T,2,Matrix{T},Tuple{UnitRange{Int},UnitRange{Int}},false}
        reduced_flag#::SubArray{Bool,2,BitMatrix,Tuple{UnitRange{Int},UnitRange{Int}},false}
    end

    """
    full_image_flag(filtered_im::FilteredImage)

Returns the BitMatrix flag of filtered pattern in the whole image.

Can be used as index matrix in the full image:

 `filtered_image.full.initial[full_image_flag(filtered_image)]` will return 
 all elements which belong to the pattern

"""
function full_image_flag(filtered_im::FilteredImage,::Type{T}=Matrix{Bool}) where T<:FlagMatrix 
        flag = T(undef,filtered_im.full.sz...)
        fill!(flag,false)
        @. flag[filtered_im.region_indices]=true
        return flag
    end   
    reduced_image_flag(fim::FilteredImage) =  copy(fim.reduced_flag)  #creates new matrix
    reduced_image(fim::FilteredImage) = copy(fim.reduced)

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
    filter_image(imag::AbstractMatrix,c::CentredObj;external=false)

Filters image according to centered object creating new image
if external  is true than as a filtering flag the inverse of centered object image is taken
"""
filter_image(imag::AbstractMatrix,c::CentredObj;external=false) = filter_image!(copy(imag),cent_to_flag(c,size(imag),external= !external))
filter_image(imag,flag::FlagMatrix) =  filter_image!(copy(imag),flag)
filter_image(imag::RescaledImage;label::Int = 1) = filter_image(imag,marker_image(imag),label=label)
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

In-place filtering of [`RescaledImage`](@ref), filtered object is wrapped around the input RescaledImage
"""
    function filter_image!(imag::RescaledImage,
        external_region_flag::FlagMatrix)
        
        external_region = view(imag.initial,external_region_flag)
        @. external_region=0.0
        @. external_region_flag = !external_region_flag
        region_area_indices = findall(external_region_flag)
        min_ind,max_ind = extrema(region_area_indices)
        @. imag.im = imag.initial
        rescale!(imag.im,inverse_intensity=imag.inverted)
        square_view = @view  imag.initial[min_ind[1]:max_ind[1],min_ind[2]:max_ind[2]]
        square_view_flag =@view external_region_flag[min_ind[1]:max_ind[1],min_ind[2]:max_ind[2]]
        return FilteredImage(imag,
                region_area_indices,
                square_view, 
                square_view_flag)       
    end
    filter_image!(imag::RescaledImage,c::CentredObj;external::Bool=false) = filter_image!(imag,cent_to_flag(c,imag.sz,external= !external))
    filtered_mean(fltrd::FilteredImage) = Statistics.mean(view(fltrd.reduced[fltrd.reduced_flag]))
    filtered_std(fltrd::FilteredImage) = Statistics.std(view(fltrd.reduced[fltrd.reduced_flag]))
