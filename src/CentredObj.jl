# this code should not be dependent on custom image types
export CentredObj,CircleObj,SquareObj,RectangleObj,
c_view,line_within_mask,along_line_distribution,
within_mask_line_points_distribution,mean_within_mask,
along_mask_line_distribution,radial_distribution
"""
    `CentredObj` is a sort of region of interest (ROI) marker object with coordinates and dimensions. 

`CentredObj` has centre coordinates and dimensions, object's center can be 
anywhere with respect to the image indices. ROI also has one or more size parameters (in pixels)
coordinates of centre are equal to CartesianIndices, first element is the row index, 
the second element is the column index!! (y and x coordinate)
This is opposite to the ImageDraw, where first Point coordinate corresponds to the column index and 
the second one to the row index. `CentredObj` can also be used as for indexing image[c] - returns all
elements of image within c, `image[c]=x` sets all elements of image to the values of x, x firstindex 
should be 1. `CentredObj` can also be used to set all image points within the ROI to a single value.
e.g. `image[c] = 30` 

By default we have three `CentredObj` types in package: 

[`CircleObj`](@ref), [`SquareObj`](@ref) and [`RectangleObj`](@ref)

To impement `CentredObj` abstraction one needs to implement:

[`side`](@ref) - returns what is assumed to be the single number characterisation of the dimensions (by default the maximum value of dimensions(c) is taken)

[`area`](@ref) - returns the CentredObj area, used in statistics calculation and in default two objs compirason see [`isless(c1::CentredObj,c2::CentredObj)`](@ref)

[`perimeter`](@ref) - Object's perimeter

[`Base.size`](@ref) - should always return the minimum size of frame which wrappes the centred obj, e.g. for circle it is [diameter,diemater]

[`is_within`](@ref) - function to check if inds are within the `CentredObj`

[`line_within_mask`](@ref) - function to check if all line points are within the `CentredObj`

[`fill_x0!`](@ref) - function to fill the optimization starting vector during `CentredObj` 
fitting the image

[`convert_to_drawable`](@ref) fucntion to convert the [`CentredObj`](@ref) to a drawable obj for `ImageDraw`

[`parnumber`](@ref) function which returns the number of parameters     

"""
abstract type CentredObj end 

"""
    parnumber(::Type{T}) where T<:CentredObj

Returns total number of parameters needed to create new object
"""
function parnumber(::Type{T}) where T<:CentredObj    error(DomainError(T,"Undefined parnumber method")) end

"""
    (::Type{T})() where T<:CentredObj

Empty object constructor
"""
function (::Type{T})() where T<:CentredObj
    obj_from_vect(T,Vector{Int}(undef,parnumber(T)))
end
"""
    copyobj(c::T) where T<:CentreObj

Copies the [`CentredObj`](@ref) creating new instance
"""
function copyobj(c::T) where T<:CentredObj 
    return obj_from_vect(T,[c.center...,c.dimensions...])
end    
"""
    Base.length(c::CentredObj)

Total number of values needed to create [`CentredObj`](@ref) of specified type
"""
function Base.length(::T) where T<:CentredObj
    return parnumber(T)
end 

function Base.show(io::IO, c::C) where C<:CentredObj
    print(io, "$(name(c)) with center at ($(center(c)[1]) , $(center(c)[2])) and size  ($(join(string.(dimensions(c)),",")))")
end
"""
    is_within(c::CentredObj,_)

Function to check if indices are within [`CentredObj`](@ref)
"""
is_within(c::CentredObj,_)=error(DomainError(typeof(c),"no `is_within` implementation"))
dimensions(c::CentredObj)=  getfield(c,:dimensions)
"""
    side(c::CentredObj)

Returns the side of `CentredObj`
"""
side(c::CentredObj) = maximum(dimensions(c))
name(c::CentredObj) =  "No name"
"""
    is_within(c::CentredObj,i::CartesianIndex)

`CartesianIndex` support
"""
function is_within(c::CentredObj,i::CartesianIndex) 
    return is_within(c,SVector(Tuple.(i)))
end 
"""
    center(c::CentredObj)

Vector of object center coordinates (MVector)
"""
center(c::CentredObj)=  c.center

"""
    revcentre(c::CentredObj)

Object coordinates in reverse order
"""
revcentre(c::CentredObj) = reverse(center(c))
resize_proportional!(c::CentredObj;percents) = begin
    @assert 0<percents<=100 "Size multiplier should be from 0 to 100"
    c.dimensions .*= precents/100
end
"""
    is_within_iterator(img::AbstractMatrix,c::CentredObj)

Iterator over all `CartesianIndices` within the `img` which are within the CentredObj `c`
"""
function is_within_iterator(img::AbstractMatrix,c::CentredObj)
    return Iterators.filter(i->is_within(c,i),keys(img))
end

"""
    c_view(img::AbstractMatrix,c::CentredObj)

Returns the view of `img` elements which are within the `CentredObj`
"""
function c_view(img::AbstractMatrix,c::CentredObj)
    return @view img[cent_to_flag(c,size(img))]
end
"""
    Base.getindex(img::AbstractMatrix,c::CentredObj)

`CentredObj` can be used for matrix indexing, `image[centred_object]` - returns the vector 
of temperatures of all points of image lying within the `centred_object` of `CentredObj`
"""
function Base.getindex(img::AbstractMatrix,c::CentredObj)
    return map(i->getindex(img,i),is_within_iterator(img,c))
end
"""
    Base.setindex!(img::Matrix,x::Array,c::CentredObj)

img[c]=x assignes all x elements to the elements of `img` with indices lying within the `CentredObj` c
"""
function Base.setindex!(img::Matrix,x::Array,c::CentredObj)
    @assert firstindex(x)==1 # check if it is not offset array
    for (ix,iimg) in enumerate(is_within_iterator(img,c))
        img[iimg] = x[ix]
    end
    return nothing
end
"""
    Base.setindex!(img::Matrix{T},x::Number,c::CentredObj) where T

Setting all elements within the `CentredObj` to a single value
"""
function Base.setindex!(img::Matrix{T},x::Number,c::CentredObj) where T
    x_T = T(x) 
    for i in is_within_iterator(img,c)
        img[i] = x_T
    end
    return nothing
end

"""
    shift!(c::CentredObj,ind::CartesianIndex{2})

Relative shift of centred object center
"""
function shift!(c::T,ind::CartesianIndex{2}) where T<:CentredObj
    c.center[1] += ind[1]
    c.center[2] += ind[2]
    return nothing
end
"""
    cent_to_flag(c::CentredObj,sz::Tuple{Int,Int};external=false)

Converts CentredObj to bitmatrix  of size sz
"""
function cent_to_flag(c::T,sz::Tuple{Int,Int};external=false) where T<:CentredObj
    external_part_flag = BitMatrix(undef,sz...)
    return external ? fill_im_external!(external_part_flag,c) : fill_im!(external_part_flag,c)
end

"""
    cent_to_flag(::Type{T},c::CentredObj,sz::Tuple{Int,Int};external=false) where T<:FlagMatrix

Converts centred obj to BitMatrix or the Matrix of bool see `FlagMatrix
"""
function cent_to_flag(::Type{T},c::CentredObj,sz::Tuple{Int,Int};external=false) where T<:FlagMatrix
    external_part_flag = T(undef,sz...)
    return external ? fill_im_external!(external_part_flag,c) : fill_im!(external_part_flag,c)
end
"""
    line_within_mask(c::CentredObj,ang::Float64,line_length::Int)

Function returns endpoint of the line lying fully within the mask  - tuple of four point which can be 
directly splatted to the along_line_distribution

ang - angle in degrees 

line_length - the length of line   
"""
function line_within_mask(c::CentredObj,::Float64,::Int)  DomainError(typeof(c),"no implementation") end
"""
    area(c::CentredObj)

Ealuates the surface area in pixels
"""
area(c::CentredObj)= error(DomainError(typeof(c),"no `area` implementation"))
"""
    perimeter(c::CentredObj)

Returns perimeter of the object
"""
perimeter(c::CentredObj)= error(DomainError(typeof(c),"no `perimeter` implementation"))
"""
    fill_x0!(x0,im_bin::AbstractMatrix,c::CentredObj)

Fills the optimization starting vector by seraching the centre of the image `im_bin`
"""
function fill_x0!(x0,im_bin::FlagMatrix,c::CentredObj) DomainError(typeof(c),"no implementation") end
"""
    convert_to_drawable(::CentredObj)

Converts CentredObj to a drawable structure appropriate to the `ImageDraw`
draw function, polygon,ellipse see [`ImageDraw.draw`] function 
"""
function convert_to_drawable(::CentredObj) end
"""
    fill_im!(img,c::CentredObj)

Fills bitmatrix `img` in a way that all pixels which are 
within the `CentredObj` are true and false otherwise.  
"""
function fill_im!(img,c::CentredObj)
    for i in keys(img)
        inds = [k for k in Tuple.(i)]
        img[i] = is_within(c,inds)
    end
    return img
end
"""
    fill_vect!(x::AbstractVector, c::CentredObj)

Converts `CentredObj` to vector
"""
function fill_vect!(x::AbstractVector, c::CentredObj)
    x[1] = c.center[1];x[2]=c.center[2];
    x[3:end] .= c.dimensions
    return x
end

"""
    image_fill_discr(image::AbstractMatrix,c::CentredObj)

Function returns the function to evaluate the discrepancy  between 
`CentredObj` and the matrix, this function is used during the fitting procedure 
"""    
function image_fill_discr(image::AbstractMatrix,c::CentredObj)
     im_copy = copy(image)   
     return x-> begin
                image_discr(image, fill_im!(im_copy,fill_from_vect!(c,x)))
        end
end
struct ImCentDiscr{T} where T<:CentredObj
    image_target::AbstractMatrix
    image_fillable::AbstractMatrix
    sz # image size
    c::T
end
(icd::ImCentDiscr)(x::AbstractVector) = begin
    fill_from_vect!(icd.c,x)
    fill_im!(icd.image_fillable,icd.c)
    icd.image
end
# arithmetic operations on CentredObj
"""
    Base.isless(c1::CentredObj,c2::CentredObj)

By default centred objects are compared by their areas, thus the vector of centred objects can be sorted
"""
Base.isless(c1::CentredObj,c2::CentredObj) = Base.isless(area(c1),area(c2))
Base.:*(c::CentredObj,a::Number) = begin 
        c_copy = copyobj(c)
        @. c_copy.dimensions = int_floor_abs(a*c_copy.dimensions)
    return c_copy
end
Base.:*(a::Number,c::CentredObj) = c*a
Base.:/(c::CentredObj,a::Number) = begin 
            c_copy = copyobj(c)
            @. c_copy.dimensions = int_floor_abs(c_copy.dimensions/a)
    return c_copy
end
"""
    obj_from_vect(::Type{CentredObj},v::AbstractVector)

Creates object from parameters vector, first two arguments are center
point other are dimensions [center[1],center[2],dimensions[1],...]
"""
function obj_from_vect(::Type{T},v::AbstractVector) where T<:CentredObj
        c = T() # empty constructor calling
        fill_from_vect!(c, v)
        return c
end

"""
    fill_from_vect!(c::CentredObj, v::AbstractVector)

Fills CentreObj parameters from the vector [center_index_1,center_index_2,dimension_1,dimension_2,...]
"""
function fill_from_vect!(c::CentredObj, v::AbstractVector)
    @assert length(c) == length(v)
    l_d = length(c.dimensions)
    c.center[1] = int_floor(v[1])
    c.center[2] = int_floor(v[2])
    map!(int_floor_abs,c.dimensions,v[3:2+l_d])
    return c
end
#-------------------------CIRCLE-OBJ---------------------------
"""
    Circle object with defined diameter
"""
mutable struct CircleObj <:CentredObj
    center::MVector{2,Int} # central point location (indices)
    dimensions::MVector{1,Int} # side length{1,Int} # diameter
    CircleObj(center,diameter::Number) = begin
        new(MVector{2}(map(int_floor_abs, center)),MVector{1}(int_floor_abs(diameter)))
    end
    CircleObj() = new(MVector{2}(1,1),MVector{1}(1))
end

name(::CircleObj) = "Circle"
parnumber(::Type{CircleObj}) = 3
diameter(c::CircleObj) = dimensions(c)[]
perimeter(c::CircleObj) = 2*π*radius(c)
radius(c::CircleObj) = diameter(c)/2
side(c::CircleObj) = diameter(c)
area(c::CircleObj) = π*radius(c)^2
is_within(c::CircleObj,inds::AbstractVector) = begin 
    #sqrt(sum(abs2  , c.center .- inds))
    return sqrt(^(c.center[1] - inds[1],2.0) + ^(c.center[2] - inds[2],2.0)) < radius(c)
end
Base.size(c::CircleObj) = (diameter(c),diameter(c))
"""
    fill_x0!(x0,im_bin::FlagMatrix,::CircleObj)

Fills starting vector for the optimization of `CentredObj`
"""
function fill_x0!(x0,im_bin::FlagMatrix,::CircleObj)

        min_ind = findfirst(im_bin)#cartesian indices
        max_ind = findlast(im_bin)

        starting_diameter= sqrt(sum(abs2, Tuple.(max_ind - min_ind)))
        x0[end] = starting_diameter
        x0[1] =0.5*(min_ind[1] + max_ind[1])
        x0[2] = 0.5*(min_ind[2] + max_ind[2])
        #x0 .= [collect(x/2 for x in Tuple.(max_ind + min_ind))..., starting_diameter]
end    
"""
    line_within_mask(c::CircleObj,ang,line_length)

Returns two endpoints of the line lying totally inside the `CentredObj`
"""
function line_within_mask(c::CircleObj,ang,line_length) 
    ang %= 360
    line_length = line_length>diameter(c) ? radius(c) : line_length/2 
    
    lsin = int_floor(line_length*sind(ang))
    lcos = int_floor(line_length*cosd(ang))

    return        [c.center[1] -  lcos,
                   c.center[2] -  lsin,
                   c.center[1] +  lcos,
                   c.center[2] +  lsin]
end
function convert_to_drawable(c::CircleObj;fill=false,thickness::Int=-1)
    thickness_threshold = int_floor(0.35*radius(c))
    if thickness==-1|| thickness >= radius(c) || thickness<thickness_threshold
        thickness = thickness_threshold
    end
    return ImageDraw.CirclePointRadius(c.center[2],c.center[1],radius(c),thickness=thickness,fill=fill)
end
#-----------------------SQUARE-OBJ---------------------------
"""
    Square with defined center and side
"""
mutable struct SquareObj <:CentredObj
    center::MVector{2,Int}
    dimensions::MVector{1,Int} # side length
    
    SquareObj(center,side) = begin
        new(MVector{2}(map(int_floor_abs, center)),MVector(int_floor_abs(side)))
    end
    SquareObj() = new(MVector{2}(1,1),MVector{1}(1))
end

name(::SquareObj) = "Square"
parnumber(::Type{SquareObj}) = 3
area(c::SquareObj)=^(dimensions(c)[],2)
side(c::SquareObj) = dimensions(c)[]
perimeter(c::SquareObj) = 4*side(c)
Base.size(c::SquareObj) = (side(c),side(c))

is_within(c::SquareObj,inds::AbstractVector) = begin
    a = side(c)/2
    abs(c.center[1] - inds[1]) <= a && abs(c.center[2] - inds[2]) <= a
    #(c.center[1]-a) <=inds[1] <=c.center[1]+a   && c.center[2]-a <=inds[2]<=c.center[2]+a
end

function fill_x0!(x0,im_bin::FlagMatrix,::SquareObj)
    
    min_ind = findfirst(im_bin)
    max_ind = findlast(im_bin)
    x0[end] = sqrt(sum(abs2, Tuple.(max_ind - min_ind)))
    x0[1] =0.5*(min_ind[1] + max_ind[1])
    x0[2] = 0.5*(min_ind[2] + max_ind[2])

end  

function line_within_mask(c::SquareObj,ang,line_length) 
    a = side(c)
    ang %=360
    if (45<=ang<=135) || (225<=ang<=315)
        a_l = abs(a/sind(ang))
        if line_length>a_l
            line_length = a_l
        end   
        lsin_f = line_length*sind(ang)/2
        lcos_f = line_length*cosd(ang)/2
    else
        a_l = abs(a/cosd(ang))
        if line_length>a_l
            line_length = a_l
        end   
        lsin_f = line_length*sind(ang)/2
        lcos_f = line_length*cosd(ang)/2         
    end

    
    lsin = int_floor(lsin_f)
    lcos = int_floor(lcos_f)

    return      [  c.center[1]- lcos,
                   c.center[2]- lsin,
                   c.center[1]+ lcos,
                   c.center[2]+ lsin ]
end
"""
    diagonal_points(c::Union{SquareObj,CircleObj})

Returns diagonal points in row-column coordinates
"""
diagonal_points(c::Union{SquareObj,CircleObj}) = begin
    a = side(c)
    a=int_floor_fld(a,2)
    return (c.center[1]-a , c.center[2]-a, c.center[1]+a , c.center[2]+a)
end
"""
    rearranged_diagonal(c::Union{SquareObj,CircleObj})

Returns diagonal points in x-y coordinates
"""
rearranged_diagonal(c::Union{SquareObj,CircleObj}) = begin
    a = side(c)
    a=int_floor_fld(a,2)
    return (c.center[2]-a,c.center[1]-a ,c.center[2]+a, c.center[1]+a )
end
#--------------------------RECTANGLE-OBJ-----------------------

"""
    Rectangular object with defined two sides
"""
mutable struct RectangleObj <:CentredObj

    center::MVector{2,Int}
    dimensions::MVector{2,Int} # two sized
   
    RectangleObj(center,sides) = begin
        d = MVector{2,Int}(undef)
        @. d = int_floor_abs(sides)
        new( MVector{2}(map(int_floor_abs, center)),d)
    end
    RectangleObj() = new(MVector{2}(1,1),MVector{2}(1,1))
end
name(::RectangleObj) = "Rectangle"
parnumber(::Type{RectangleObj}) = 4
area(c::RectangleObj)=*(sides(c)...)
perimeter(c::RectangleObj) = 2*sum(sides(c))
sides(c::RectangleObj) = (c.dimensions[1],c.dimensions[2])
side(c::RectangleObj) = maximum(sides(c))
Base.size(c::RectangleObj) = sides(c)
is_within(c::RectangleObj,inds::AbstractVector) = begin
    (a,b) = sides(c)
    a/=2
    b/=2
    c.center[1]-a<=inds[1]<=c.center[1]+a   && c.center[2]-b<=inds[2]<=c.center[2]+b
end
diagonal_points(c::RectangleObj) = begin
    (a,b) = sides(c)
    a=int_floor_fld(a,2) # corresponds to vertical size(row index)
    b=int_floor_fld(b,2) # horizontal side (column index)
    return (c.center[1]-a , c.center[2]-b, c.center[1]+a , c.center[2]+b)
end
rearranged_diagonal(c::RectangleObj) = begin
    (a,b) = sides(c)
    a=int_floor_fld(a,2)
    b=int_floor_fld(b,2)
    return (c.center[2]-b,c.center[1]-a , c.center[2]+b, c.center[1]+a )
end
function convert_to_drawable(c::Union{RectangleObj,SquareObj};kwargs...)
    return  ImageDraw.Polygon(ImageDraw.RectanglePoints(rearranged_diagonal(c)...))
end

function fill_x0!(x0,im_bin::FlagMatrix,::RectangleObj)
    #(min_ind,max_ind) = extrema(findall(im_bin))
    min_ind = findfirst(im_bin)
    max_ind = findlast(im_bin)

    starting_a = sqrt(sum(abs2, Tuple.(max_ind - min_ind)))
    starting_b = starting_a
    x0 .= [collect(x/2 for x in Tuple.(max_ind + min_ind))...,starting_a,starting_b]
end 
function diag_ang(c::RectangleObj)
    a,b = sides(c)
    return atand(a/b)
end
function line_within_mask(c::RectangleObj,ang,line_length) 
    b,a = sides(c) # here we interchange the sides
    ang %=360
    rect_ang = diag_ang(c)
    if ((90-rect_ang)<=ang<=(180-rect_ang)) || ((270-rect_ang)<=ang<=(360-rect_ang))
        a_l = abs(a/sind(ang))
        if line_length>a_l
            line_length = a_l
        end   
        lsin_f = line_length*sind(ang)/2
        lcos_f = line_length*cosd(ang)/2
    else
        a_l = abs(b/cosd(ang))
        if line_length>a_l
            line_length = a_l
        end   
        lsin_f = line_length*sind(ang)/2
        lcos_f = line_length*cosd(ang)/2         
    end

    
    lsin = int_floor(lsin_f)
    lcos = int_floor(lcos_f)

    return      [  c.center[1]- lcos,
                   c.center[2]- lsin,
                   c.center[1]+ lcos,
                   c.center[2]+ lsin ]
end
"""
    along_line_distribution(img::AbstractMatrix{T},x0,y0,x1,y1) where T

Function evaluates matrix values distribution along the line specified by two coordinates, 
img - input image 
returns the tuple of two vectors: coordinates and values 
see `ImageDraw.bresenham` for details of finding the points of the line 

returns 
points - vector of coordinates along the line
distrib - distribution
"""
function along_line_distribution(img::AbstractMatrix{T},x0,y0,x1,y1) where T
        dx = abs(x1 - x0)
        dy = abs(y1 - y0)

        sx = x0 < x1 ? 1 : -1
        sy = y0 < y1 ? 1 : -1;

        err = (dx > dy ? dx : -dy) / 2
        points = Vector{CartesianIndex}()
        distrib = Vector{T}()
        while true
            Base.push!(points,CartesianIndex(x0,y0))
            Base.push!(distrib,img[x0,y0])
            (x0 != x1 || y0 != y1) || break
            e2 = err
            if e2 > -dx
                err -= dy
                x0 += sx
            end
            if e2 < dy
                err += dx
                y0 += sy
            end
        end
    return (points,distrib)
end



"""
    mean_within_mask(img::AbstractMatrix,c::CentredObj)

Evaluates the average temperature of all points within the `CentredObj` marker
"""
function mean_within_mask(img::AbstractMatrix,c::CentredObj)
    flt = Iterators.filter(i->is_within(c,i),keys(img))
    return Statistics.mean(i->img[i],flt)
end


"""
    add_distrib_point!(points,distrib,point,value)

Internal fucntion to add the point to distribution
"""
add_distrib_point!(points,distrib,point,value) = begin
    Base.push!(points,point)
    Base.push!(distrib,value)       
end
"""
    along_line_distribution_xiaolin_wu(img::AbstractMatrix{T}, y0, x0, y1, x1) where T

Evaluates the value matrix content along the line with endpoint coordinates x0,y0,y1,x1,
returns indices of all points. As far as Wu's algorithm returns two adjacent points
the value is evaluated as an average of two point obtained with Wu's algorithm

see  `xiaolin_wu` function from `ImageDraw` 
"""
function along_line_distribution_xiaolin_wu(img::AbstractMatrix{T}, y0, x0, y1, x1) where T
        dx = x1 - x0
        dy = y1 - y0

        swapped=false
        if abs(dx) < abs(dy)
            x0, y0 = swap(x0, y0)
            x1, y1 = swap(x1, y1)
            dx, dy = swap(dx, dy)
            swapped=true
        end
        if x1 < x0
            x0, x1 = swap(x0, x1)
            y0, y1 = swap(y0, y1)
        end
        gradient = dy / dx

        points = Vector{CartesianIndex}()
        distrib = Vector{T}()

        xend = round(Int, x0)
        yend = y0 + gradient * (xend - x0)
        xgap = rfpart(x0 + 0.5)
        xpxl0 = xend
        ypxl0 = trunc(Int, yend)
        index = swapped ? CartesianIndex(xpxl0, ypxl0) : CartesianIndex(ypxl0, xpxl0)

        #drawifinbounds!(img, index, T(rfpart(yend) * xgap))
        prev_val = img[index]
        index = swapped ? CartesianIndex(xpxl0, ypxl0 + 1) : CartesianIndex(ypxl0 + 1, xpxl0)
        #drawifinbounds!(img, index, T(fpart(yend) * xgap))
        add_distrib_point!(points,distrib,index,0.5*(img[index]+prev_val))

        intery = yend + gradient
        xend = round(Int, x1)
        yend = y1 + gradient * (xend - x1)
        xgap = fpart(x1 + 0.5)
        xpxl1 = xend
        ypxl1 = trunc(Int, yend)

        index = swapped ? CartesianIndex(xpxl1, ypxl1) : CartesianIndex(ypxl1, xpxl1)
        prev_val = img[index]

        index = swapped ? CartesianIndex(xpxl1, ypxl1 + 1) : CartesianIndex(ypxl1 + 1, xpxl1)
        add_distrib_point!(points,distrib,index,0.5*(img[index]+prev_val))


        for i in (xpxl0 + 1):(xpxl1 - 1)
            index = swapped ? CartesianIndex(i, trunc(Int, intery)) : CartesianIndex(trunc(Int, intery), i)
            prev_val = img[index]

            index = swapped ? CartesianIndex(i, trunc(Int, intery) + 1) : CartesianIndex(trunc(Int, intery) + 1, i)
            add_distrib_point!(points,distrib,index, 0.5*(img[index]+prev_val))

            intery += gradient
        end
        inds = fill(0,length(points))
        sortperm!(inds,points)
        points .=points[inds]
        distrib .=distrib[inds]
        return (points,distrib)
    end
    function swap(x, y)
        y, x
    end
    fpart(pixel::T) where {T} = pixel - T(trunc(pixel))
    rfpart(pixel::T) where {T} = oneunit(T) - fpart(pixel)
"""
    within_mask_line_points_distribution(imag::AbstractMatrix,c::CentredObj,direction_angle=0.0,line_length=10.0;use_wu::Bool=false)

Evaluates the distribution of values in `imag` matrix along the line with length `line_length` in pixels
oriented with the angle `direction_angle` in degrees  with respect to the posistive direction of oX (column index increase), 
this line lies within the roi (`CentreObj`) and goes through its center.

Function returns:

points  - vector of `CartesianIndex` of image's points lying on the line

distrib - distribution of values

line_points - endpoints of line the Tupple of (left_x,left_Y,right_x,right_y)

"""
function within_mask_line_points_distribution(imag::AbstractMatrix,c::CentredObj,
                        direction_angle=0.0,line_length=10.0;
                        use_wu::Bool=false)

    line_points = line_within_mask(c,direction_angle,line_length) 
    points_within_line!(imag,line_points)
    if use_wu
        (points,distrib) = along_line_distribution_xiaolin_wu(imag,line_points...)
    else
        (points,distrib) = along_line_distribution(imag,line_points...)
    end
    return (points,distrib,line_points)
end
"""
    along_mask_line_distribution(imag::AbstractMatrix,c::CentredObj,direction_angle=0.0,line_length=10.0;
                                                                                            length_per_pixel=1.0,
                                                                                            use_wu::Bool=false)


The same as `within_mask_line_points_distribution` but returns the line length along the coordinates within the 
image.

`line_length` - the length of line in the same units as `length_per_pixel`.
The calibration `using mm_per_pixel``, returns calibrated length along the line 

"""
function along_mask_line_distribution(imag::AbstractMatrix,c::CentredObj,
                                    direction_angle=0.0,  line_length=10.0;
                                    length_per_pixel=1.0, use_wu::Bool=false)
    
    line_length = line_length/length_per_pixel         # converting line_length to pixels                                                                            
    (points,distrib,line_points) = within_mask_line_points_distribution(imag,c,direction_angle,line_length,
                                                                                                    use_wu=use_wu)
    along_line_length = line_points_to_along_length(points,line_points)*length_per_pixel
    return (along_line_length,distrib,line_points)
end
"""
    radial_distribution(imag::AbstractMatrix,c::CentredObj,angles_range::AbstractRange,line_length;mm_per_pixel=1.0)

Calls `along_mask_line_distribution` on lines oriented with some angles range and puts the resulting 
distribution into one matrix j'th column of this matrix corresponds to the distribution along the 
line oriented with ang[j] angle. The length of line is `line_length`, if it is less than 0.0 or greater 
than the smallest dimension of `c`, than is is set to the smallest dimentsion of `c`

"""
function radial_distribution(imag::AbstractMatrix,c::CentredObj,
                            angles_range::AbstractRange;line_length=0.0,
                                            length_per_pixel=1.0,
                                            use_wu::Bool=false)
    # first calling to obtain the length
    line_length = line_length <=0.0 ? minimum(c.dimensions) : line_length 
    (first_columnn_along_line,_,) = along_mask_line_distribution(imag,c,0.0, line_length;length_per_pixel=length_per_pixel,use_wu=use_wu)
    # 
    points_number = length(first_columnn_along_line)
    angles_number = length(angles_range)
    radial_distrib_matrix = fill(NaN,points_number,angles_number)#Matrix{Float64}(undef,points_number,angles_number)

    extrapolation_bc = Interpolations.Line()
    #extrapolation_bc = Interpolations.Flat()
    Threads.@sync for (i,α) in enumerate(angles_range) # circle over line rotation angle 
            Threads.@spawn begin 
                (along_line_length,distrib,) = along_mask_line_distribution(imag,c,α, line_length;length_per_pixel=length_per_pixel,use_wu=use_wu)
                w_distr = @view radial_distrib_matrix[:,i]
                along_line_distrib = LinearInterpolation(along_line_length,distrib,extrapolation_bc = extrapolation_bc)(first_columnn_along_line) 
                copyto!(w_distr,along_line_distrib)
            end
    end
    return (first_columnn_along_line,radial_distrib_matrix)
end