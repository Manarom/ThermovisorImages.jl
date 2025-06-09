
"""
Special type to work with segmentated image 
"""
struct MarkeredImage
    markers::Matrix{Int} # this stores the matrix of markers
    ViewsVect::Vector{Vector{CartesianIndex{2}}} # stores a vector of indices per each pattern
    MarkeredImage(markers::Matrix{Int}) = begin
        flag = similar(markers,Bool)
        max_label = maximum(markers)
        if max_label<1
            return new(markers,Vector{Vector{CartesianIndex{2}}}([])) 
        end
        map!(i->i==1 ,flag, markers)
        ViewsVect = Vector{Vector{CartesianIndex{2}}}(undef,max_label)
        ViewsVect[1] = findall(flag)
        for i in 2:max_label
            map!(m->m==i ,flag, markers)
            ViewsVect[i] = findall(flag)
        end
        return new(markers,ViewsVect)    
    end
end
"""
Base.length(m::MarkeredImage)

Return the total number of patterns in the markerred image [`MarkeredImage`](@ref)
"""
Base.length(m::MarkeredImage) = length(m.ViewsVect)

"""
Base.size(m::MarkeredImage)

The image size in pixels
"""
Base.size(m::MarkeredImage) = size(m.markers)

#indexing markered image
Base.lastindex(m::MarkeredImage) = length(m)
Base.getindex(m::MarkeredImage,i) = m.markers[m.ViewsVect[i]]
Base.setindex!(m::MarkeredImage,val::Int,i::Int) = fill!(m_view(m,i),val)
m_view(m::MarkeredImage,i::Int) = view(m.markers,m.ViewsVect[i])


"""
areas(m::MarkeredImage)

Vector of patterns areas (number of pixels within the pattern)
"""
areas(m::MarkeredImage) = map(length, m.ViewsVect)

pattern_diagonal(m::MarkeredImage,i::Int) = extrema(m.ViewsVect[i])
"""
flag(m::MarkeredImage,i::Int)


"""
function flag(m::MarkeredImage,i::Int)
    fl = BitMatrix(undef,size(m.markers))
    return flag!(fl,m,i)
end
"""
external_flag(m::MarkeredImage,i::Int)

Inversed version of [`flag`](@ref)
"""
function external_flag(m::MarkeredImage,i::Int)
    fl = BitMatrix(undef,size(m.markers))
    return flag!(fl,m,i,negate=true)
end
"""
flag!(fl::FlagMatrix,m::MarkeredImage,i::Int;negate::Bool=false)

Fills the fl matrix (Bitmatrix or Matrix{Bool}) with the same size 
as the entire image with all elements set to zero, except the pixels of `i`'th pattern
"""
function flag!(fl::FlagMatrix,m::MarkeredImage,i::Int;negate::Bool=false)
        fill!(fl,negate)
        f = @view fl[m.ViewsVect[i]]
        fill!(f,!negate)
        return fl
end
"""
reduced_flag(m::MarkeredImage,i::Int)

returns flag matrix shrinked to the size of the pattern
"""
function shrinked_flag(m::MarkeredImage,i::Int)
        (min_i, max_i) = pattern_diagonal(m,i) #CartesianIndices of maximum and minimum indices of the pattern
        #sz = 
        dinds = max_i - min_i + CartesianIndex(1,1)#indices difference
        fl = falses(Tuple.(dinds))# creating the size of obj
        for c in m.ViewsVect[i]
            i_in = c - min_i + CartesianIndex(1,1)
            fl[i_in] = true#m.markers[c]
        end
        return (fl,min_i,max_i)
end

"""
sort_reduce!(m::MarkeredImage;total_number::Int = -1,descending::Bool=true)

Sorts markers by total area and reduces the total number of patterns if the maximum label
is less then `total_number` value

Input arguments:
m - [`MarkeredImage`](@ref)
total_number - number of 
"""
function sort_reduce!(m::MarkeredImage;total_number::Int = -1,descending::Bool=true)
        max_label = length(m)
        is_reduce = !(total_number <=0 || total_number >max_label) 
        if !is_reduce
            total_number = max_label
        end
        inds = Vector{Int}(undef,max_label)
        areas_array = areas(m) # areas order is the same as patterns order 
        sortperm!(inds,areas_array,rev=descending)
        permute!(m.ViewsVect,inds) # now indices vectors are in the same order as the areas 
        # first array of indices corresponds to the patter with the maximal area
        for i = 1:total_number
            m[inds[i]] = i
        end
        if is_reduce 
            for i= total_number+1:max_label
                m[inds[i]] = 0
            end
            deleteat!(m.ViewsVect,total_number+1:max_label)
        end
    return m
end
sort_markers!(m::MarkeredImage,descending::Bool=true) = sort_reduce!(m,total_number=-1,descending=descending)

"""
    count_separate_patterns(markers::Matrix{Int})

This function takes matrix of markers see [`marker_image`](@ref) and calculates the number of separate patterns
"""
function   count_separate_patterns(markers::Matrix{Int})
    return maximum(markers)
end

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


function Base.show(io::IO,m::MarkeredImage)
    print(io, "Markered image of size $(size(m)) with $(length(m)) markered patterns")
end