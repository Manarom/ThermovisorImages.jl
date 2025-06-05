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
 