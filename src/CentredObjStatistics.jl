
export CentredObjCollectionStat, std_within_mask,
radial_distribution_statistics,DistributionStatistics,
angular_distribution_statistics,
generate_random_objs

max_side(v) = findmax(side,v)
max_area(v) = findmax(area,v)
max_perimeter(v) = findmax(perimeter,v)

min_side(v) = findmin(side,v)
min_area(v) = findmin(area,v)
min_perimeter(v) = findmin(perimeter,v)

mean_side(v) = Statistics.mean(side,v)
mean_area(v) = Statistics.mean(area,v)
mean_perimeter(v) = Statistics.mean(perimeter,v)

std_side(v) = Statistics.std(map(side,v))
std_area(v) = Statistics.std(map(area,v))
std_perimeter(v) = Statistics.std(map(perimeter,v))

hist_side(v,nbins::Int=-1) = nbins<=0 ? fit(Histogram,map(side,v)) : fit(Histogram,map(side,v),nbins=nbins)
hist_area(v,nbins::Int=-1) = nbins<=0 ? fit(Histogram,map(area,v)) : fit(Histogram,map(area,v),nbins=nbins)
hist_perimeter(v,nbins::Int=-1) = nbins<=0 ? fit(Histogram,map(perimeter,v)) : fit(Histogram,map(perimeter,v),nbins=nbins)
"""
    CentredObjCollectionStat

Calculates the following statistics on `CentredObj`s collection:

minimal, maximal, mean and standard deviation of sides,areas and perimeters of all `CentredObj`
objects in the collection `v`.
"""
struct CentredObjCollectionStat{T}
    N::Int
    maxs::T
    mins::T
    means::T
    stds::T
    hists::T
    """
    CentredObjCollectionStat(v::Vector{C};nbins::Int=-1) where C<:CentredObj

Calculates statistics on `CentredObj`s vector:

minimal, maximal, mean and standard deviation of sides,areas and perimeters of all `CentredObj`
objects in the collection `v`.
"""
    CentredObjCollectionStat(v::Vector{C};nbins::Int=-1) where C<:CentredObj= begin
            new{NamedTuple{(:side,:area,:perimeter)}}(length(v),
                (max_side(v),max_area(v),max_perimeter(v)),
                (min_side(v),min_area(v),min_perimeter(v)),
                (mean_side(v),mean_area(v),mean_perimeter(v)),
                (std_side(v),std_area(v),std_perimeter(v)),
                (hist_side(v,nbins),hist_area(v,nbins),hist_perimeter(v,nbins))
            )
    end
end

function Base.show(io::IO,s::CentredObjCollectionStat)
    println(io,)
    println(io,"""Statistics summary on $(s.N) rois:
                side = $(s.means.side) ±  $(s.stds.side) in range from $(s.mins.side[1]) to  $(s.maxs.side[1])
                area = $(s.means.area) ±  $(s.stds.area) in range from  $(s.mins.area[1]) to  $(s.maxs.area[1])
                perimeter = $(s.means.perimeter) ±  $(s.stds.perimeter) in range from  $(s.mins.perimeter[1]) to  $(s.maxs.perimeter[1])
    """)
end

"""
    std_within_mask(img::AbstractMatrix, c::CentredObj)

Evaluates standard deviation of temperature for all points within the `CentredObj` marker
"""
function std_within_mask(img::AbstractMatrix, c::CentredObj)
    flt = Iterators.filter(i->is_within(c,i),keys(img))
    mpper = Iterators.map(i->img[i],flt)
    return Statistics.std(mpper)
end
"""
    radial_distribution_statistics(along_length_coordinate::AbstractVector,
        distrib::AbstractVecOrMat;
        max_length=-1.0,min_length=-1.0)

Fills [`DistributionStatistics`](@ref) object for radial distribution

Optional:

`max_length` - maximal value of along_length_coordinate to be includet in to the statistics evaluation

`min_length` - minimal value of along_length_coordinate to be includet in to the statistics evaluation

"""
function radial_distribution_statistics(along_length_coordinate::AbstractVector,
        distrib::AbstractVecOrMat;
        max_length=-1.0,min_length=-1.0)
            @assert(length(along_length_coordinate)==size(distrib,1),"Vector of coordinate should have the same length as the number of distrib rows")
            not_nan_flag = _inbounds_flag(along_length_coordinate,distrib,max_length,min_length)
            L = along_length_coordinate[not_nan_flag]
            D = @view distrib[not_nan_flag,:]
            return DistributionStatistics(L,D)
end
"""
    _inbounds_flag(L,D,max_length,min_length)

Unsafe version of check! number of  rows in `D` should be the same as the number of 
elements in `L`
Returns Bool flag of all row not containing NaN's and lying within the min_length to max_length range
"""
function _inbounds_flag(L,D,max_length,min_length)
    not_nan_flag = Vector{Bool}(undef,size(D,1)) 
    if max_length>0 && max_length < maximum(L)
        if  min_length>0 && min_length>minimum(L)
            for (i,r) in enumerate(eachrow(D))
                @inbounds not_nan_flag[i] = !any(isnan,r) && L[i]<=max_length && L[i]>=min_length
            end
        else
            for (i,r) in enumerate(eachrow(D))
                @inbounds not_nan_flag[i] = !any(isnan,r) && L[i]<=max_length
            end
        end
    else
        for (i,r) in enumerate(eachrow(D))
            @inbounds not_nan_flag[i] = !any(isnan,r)
        end
    end
    return not_nan_flag
end
"""
    DistributionStatistics

Basic descriptive statistics (mean and standard deviation) on temeprature distribution 
of value for a given `coordinate`,   `D ` is the matrix of distribution, where each
column corresponds to sample, thus rows number of `D` should be the same as the length 
of `coordinate`, columns number is the number of samples.

"""
struct DistributionStatistics{T}
    coordinate::T
    mean_D::T
    std_D::T
    samples_number::Int

"""
    DistributionStatistics(coordinate::AbstractVector,D::AbstractMatrix)

Basic descriptive statistics (mean and standard deviation) on temeprature distribution 
of value for a given `coordinate`,   `D ` is the matrix of distribution, where each
column corresponds to sample, thus rows number of `D` should be the same as the length 
of `coordinate`, columns number is the number of samples.

"""
DistributionStatistics(coordinate::AbstractVector,D::AbstractMatrix) = begin
        len = length(coordinate)
        len != size(D,1) ? throw(DomainError("Inappropriate size of D and coordinates")) : nothing
        mean_D = similar(coordinate)
        Statistics.mean!(mean_D,D)
        std_D = vec(Statistics.stdm(D,mean_D,dims=2))
        samples_number = size(D,2)
        return new{typeof(coordinate)}(copy(coordinate),mean_D,std_D,samples_number)   
    end     
end
 """
    eval_bounds(DS::DistributionStatistics;is_use_student::Bool=true)

Evaluates confidence bounds
"""
function eval_bounds(DS::DistributionStatistics;is_use_student::Bool=true,probability::Float64=0.95)
    t_value =  is_use_student ?  student_coefficient(DS.samples_number,probability) : 1.0
    l_b = copy(DS.mean_D)
    u_b = copy(DS.mean_D)
    @. l_b -= t_value*DS.std_D
    @. u_b += t_value*DS.std_D
    return (l_b,u_b,t_value)
 end
"""
    angular_distribution_statistics(angles,along_length_coordinate,distrib;
                            max_length=-1.0,min_length=-1.0)

Fills [`DistributionStatistics`](@ref) object for angular

Optional:

`max_length` - maximal value of along_length_coordinate to be includet in to the statistics evaluation

`min_length` - minimal value of along_length_coordinate to be includet in to the statistics evaluation
"""
function angular_distribution_statistics(angles,along_length_coordinate,distrib;
                            max_length=-1.0,min_length=-1.0)
    
    not_nan_flag = _inbounds_flag(along_length_coordinate,distrib,max_length,min_length)
    #L = @view along_length_coordinate[not_nan_flag]
    D =transpose( @view distrib[not_nan_flag,:])
    return DistributionStatistics(angles,D)
end
"""
    points_within_line!(imag::AbstractMatrix,line_points::AbstractVector)

Forces all line points to lie within the possible region according to the image size
"""
function points_within_line!(imag::AbstractMatrix,line_points::AbstractVector)
    sz = size(imag)
    for (ind,l) in enumerate(line_points)
        if l<=0
            line_points[ind] = 1
        else
            s = isodd(ind) ? sz[1] : sz[2]
            if  l>s
                line_points[ind] = s
            end
        end
    end
    return line_points
end

"""
    line_points_to_along_length(along_line_points::Vector{T},line_points) where T

Converts Cartesian indices of `along_line_points` to the length along line
"""
function line_points_to_along_length(along_line_points::Vector{T},line_points) where T
    line_start = T(Tuple(line_points[1:2]))
    length_along_line = Vector{Float64}(undef,length(along_line_points))
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


    """
    generate_random_objs(::Type{T},centers_range::NTuple{2,R},obj_number::Int,dimension_range::R) where {T<:CentredObj,R<:StepRange{Int,Int}}


Generates `obj_number` of objects of `CentredObj` subtype  specified by the first argument, 
Second argument `centers_range` is a tuple of `StepRange` s for x and y coordinates of center 
to be taken from randomly, and `dimension_range` is the range for dimensions
"""
function generate_random_objs(::Type{T},centers_range::NTuple{2,R},obj_number::Int,dimension_range::R) where {T<:CentredObj,R<:StepRange{Int,Int}}
        @assert obj_number>0 "The number of patterns should be greater than zero"
        dim_number = parnumber(T) - 2 
        rnd_centre() = [rand(centers_range[1]);rand(centers_range[2])] # random center positions
        rnd_diam() = rand(dimension_range,dim_number) # random diameter generator
        # filling the vector of initial ROIs
        return  [obj_from_vect(T,vcat(rnd_centre(),rnd_diam())) for _ in 1:obj_number]
    end

"""
    generate_random_objs!(im::Matrix{Float64},::Type{T},obj_number::Int,dimension_range::R) where {T<:CentredObj,R<:StepRange{Int,Int}}

Generates random objects on the image `im` see [`generate_random_objs`](@ref)
"""
function generate_random_objs!(im::Matrix{Float64},::Type{T},obj_number::Int,dimension_range::R) where {T<:CentredObj,R<:StepRange{Int,Int}}
    centers_range = (1:1:size(im,1),1:1:size(im,2))
    objs = generate_random_objs(T,centers_range,obj_number,dimension_range)
    (min_im,max_im) = extrema(im)
    if min_im==max_im
        min_im=0
        max_im=1
    end
    areas = map(area,objs)
    (min_a,max_a) = extrema(areas)
    if min_a==max_a
        conv=(x)->x
    else
        conv = (x)-> min_im + (max_im - min_im)*(x-min_a)/(max_a - min_a)
    end
    for (i,c) in enumerate(objs)
        im[c] = conv(areas[i])
    end
end