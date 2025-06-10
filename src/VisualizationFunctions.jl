export draw!,draw,to_rgb,draw_line_within_mask,
draw_line_within_mask!,change_default_colorscheme,
change_default_roi_color
import ImageDraw
"""
    draw!(image::Matrix{Float64},c::CentredObj;fill=false,thickness::Int=-1,
                                        roi_color::RGB{Float64}=DefRoiColor[], 
                                        color_scheme::Symbol=:none,
                                        show_cross=true,kwargs...)


Draws `CentreObj` `c` inside the `image`.

image - image

c - object 

fill - if true the interior of the object will be filled 

thickness - the thickness of the object's frame

color - frame and filling color 
"""
function draw!(image::Matrix{Float64},c::CentredObj;fill=false,thickness::Int=-1,
                                        roi_color::RGB{Float64}=DefRoiColor[], 
                                        color_scheme::Symbol=:none,
                                        show_cross=true,kwargs...) 

        rgbim = to_rgb(image,color_scheme=color_scheme)
        #im_pic = ImageDraw.draw!(rgbim,LineTwoPoints(points_inds...), RGB{Float64}(1,0,0))             
        return   draw!(rgbim,c;
            fill=fill,thickness=thickness, 
            roi_color=roi_color, show_cross=show_cross,kwargs...)   
    end
    draw(image::Matrix{Float64};color_scheme::Symbol=:none) = to_rgb(image,color_scheme=color_scheme)
    draw(image::RescaledImage;color_scheme::Symbol=:none) = draw(image.im,color_scheme=color_scheme)
    draw(image::FilteredImage;color_scheme::Symbol=:none,draw_reduced::Bool=false) = draw_reduced ? draw(reduced_image(image),color_scheme=color_scheme) : draw(image.full,color_scheme=color_scheme) 
    draw(image::RescaledImage,c::CentredObj) = begin
        im_rgb = to_rgb(image)
        draw!(im_rgb,c)
        return im_rgb
    end
    draw(image::RescaledImage,c::Vector{T}) where T<:CentredObj= begin
        im_rgb = to_rgb(image)
        for ci in c
            draw!(im_rgb,ci)
        end
        return im_rgb
    end
    function draw!(rgbim::Matrix{RGB{Float64}},
                    c::CentredObj;
                    fill=false,
                    thickness::Int=55,
                    roi_color::RGB{Float64}=DefRoiColor[],
                    show_cross::Bool=true,
                    kwargs...) 

        ImageDraw.draw!(rgbim,convert_to_drawable(c,fill=fill,thickness=thickness), roi_color; kwargs...)
        
        if show_cross 
            cross_length = int_floor(minimum(side(c))/3)
            if cross_length <5
                cross_length = side(c)
            end
            ImageDraw.draw!(rgbim, ImageDraw.Cross(ImageDraw.Point(revcentre(c)...), cross_length), roi_color)
        end    
        return rgbim
    end

    """
    to_rgb(image::Matrix{Float64};color_scheme::String="")

Converts matrix to rgb - martix by applyting the color scheme 
from `ColorSchemes` package. 
"""
function to_rgb(image::Matrix{Float64};
                color_scheme::Symbol=:none)

        return to_rgb(RescaledImage(image),color_scheme=color_scheme)
        #collect(colorview(RGB, permuteddimsview(rgbimg_3D,(3,1,2))) )
    end
function to_rgb(image::RescaledImage{Float64},color_scheme::Symbol=:none) 
    haskey(colorschemes,color_scheme) ? change_default_colorscheme(color_scheme) : nothing
    return get.(DefColorScheme,image.im)
end
    """
    draw(c::CentredObj;kwargs...)

Returns `CentredObj` image of minimal possible size
"""
function draw(c::CentredObj;kwargs...) 
        (x_left,y_left,x_right,y_right) = abs.(diagonal_points(c))
        image = fill(0.0, [y_right + y_left, x_right + x_left]...)
        image[1,1] = 0.001            
        return draw!(image,c;kwargs...)
end
"""
    draw_line_within_mask(image::Matrix{Float64},c::CentredObj,
                    ang,length;thickness::Int=55,
                    roi_color::RGB{Float64}=DefRoiColor[], 
                    color_scheme::Symbol=:none,kwargs...)

Draws straight line on the `image` converted to RGB, all points are
located within the `CentredObj`. This line goes through the centre of `c`
and oriented with the angle `ang` in degrees with positive direction  - counterclockwise.

"""
function draw_line_within_mask(image::Matrix{Float64},c::CentredObj,
                    ang,length;thickness::Int=55,
                    roi_color::RGB{Float64}=DefRoiColor[], 
                    color_scheme::Symbol=:none,kwargs...) 
    return draw_line_within_mask!(to_rgb(image,color_scheme=color_scheme),c,ang,length;thickness=thickness,
                                roi_color=roi_color, kwargs...)
end
function draw_line_within_mask!(rgbim::Matrix{T},
                                c::CentredObj,ang,length;
                                thickness::Int=55,
                                roi_color::T=T(0,1,0), kwargs...) where T<:RGB{Float64}

        line_coords = line_within_mask(c,ang,length)
        #we should interchange the order of coordinates according to the ImageDraw demands
        ImageDraw.draw!(rgbim,ImageDraw.LineSegment(line_coords[2],line_coords[1],line_coords[4],line_coords[3]), roi_color)
        return rgbim
end
 
"""
    change_default_colorscheme(new_scheme::Symbol)

Change colorschee which is used by default to convert matrices to rgb images
"""
function change_default_colorscheme(new_scheme::Symbol)
    @assert haskey(colorschemes,new_scheme) "This colorscheme is not supported see keys(ColorSchemes.colorscheme)"
    DefColorScheme[] = colorschemes[new_scheme]
end

"""
    change_default_roi_color(color::RGB{Float64})

Changes default color to visualize the `CentredObj`
"""
function change_default_roi_color(color::RGB{Float64})
    DefRoiColor[] = color
end