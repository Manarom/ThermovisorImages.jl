export draw!,draw,to_rgb,draw_line_within_mask,
draw_line_within_mask!,change_default_colorscheme,
change_default_roi_color,radial_distribution_statistics_plot,angular_distribution_statistics_plot,along_line_distribution_plot

"""
    draw(image::Matrix{Float64},c::CentredObj;fill=false,thickness::Int=-1,
                                        roi_color::RGB{Float64}=DefRoiColor[], 
                                        color_scheme::Symbol=:none,
                                        show_cross=true,kwargs...)


Converts `image` to RGB (if not all values of the `image` are within [0,1]  - rescales) \n
And draws `CentreObj` `c` inside the `image` (does not affect the values of `image`)\n
For `kwargs` see [`draw!`](@ref)
"""
function draw(image::Matrix{Float64},c::CentredObj;color_scheme::Symbol=:none,kwargs...) 

        rgbim = to_rgb(image,color_scheme=color_scheme)
        #im_pic = ImageDraw.draw!(rgbim,LineTwoPoints(points_inds...), RGB{Float64}(1,0,0))             
        return   draw!(rgbim,c;kwargs...)   
    end
    draw(image::Matrix{Float64};color_scheme::Symbol=:none) = to_rgb(image,color_scheme=color_scheme)
    draw(image::RescaledImage;color_scheme::Symbol=:none) = to_rgb(image,color_scheme=color_scheme)
    draw(image::FilteredImage;color_scheme::Symbol=:none,draw_reduced::Bool=false) = draw_reduced ? draw(reduced_image(image),color_scheme=color_scheme) : draw(image.full,color_scheme=color_scheme) 
    """
    draw(image::RescaledImage,c::CentredObj)

Converts image to RGB and draws `CentredObj` in it.
For `kwargs` see [`draw!`](@ref)
"""
draw(image::RescaledImage,c::CentredObj;kwargs...) = begin
        im_rgb = to_rgb(image)
        draw!(im_rgb,c;kwargs...)
        return im_rgb
    end
    """
    draw(image::RescaledImage,c::Vector{T};kwargs...) where T<:CentredObj

Converts `image` to RGB and draws multiple objects
"""
draw(image::RescaledImage,c::Vector{T};kwargs...) where T<:CentredObj= begin
        im_rgb = to_rgb(image)
        for ci in c
            draw!(im_rgb,ci;kwargs...)
        end
        return im_rgb
    end
    """
    draw!(rgbim::Matrix{RGB{Float64}},
                    c::CentredObj;
                    fill=false,
                    thickness::Int=55,
                    roi_color::RGB{Float64}=DefRoiColor[],
                    show_cross::Bool=true,
                    kwargs...)

Draws `c` inside the image `rgbim` (modified) \n
Keyword arguments: \n
fill - if true fills the object \n
thickness - thickness of `c` boundary \n
(this two arguments are transfered to [`convert_to_drawable`](@ref)) \n
roi_color - RGB color of ROI boundary and filling  \n
show_cross - if true the cross is displayed \n
kwargs - Other keyword arguments are transfered directly
to the `ImageDraw.draw!` function
"""
function draw!(rgbim::Matrix{RGB{Float64}},
                    c::CentredObj;
                    fill=false,
                    thickness::Int=-1,
                    roi_color::RGB{Float64}=DefRoiColor[],
                    show_cross::Bool=true,
                    kwargs...) 

        ImageDraw.draw!(rgbim,convert_to_drawable(c,fill=fill,thickness=thickness), roi_color; kwargs...)
        
        if show_cross 
            cross_length = int_floor(side(c)/3)
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
        (min,max) = extrema(image)
        is_scaled =  min >=0.0 && max <=1 #if it is already scaled 
        return   is_scaled ? _to_rgb(image,color_scheme=color_scheme) : to_rgb(RescaledImage(image),color_scheme=color_scheme)
        #collect(colorview(RGB, permuteddimsview(rgbimg_3D,(3,1,2))) )
    end
"""
    _to_rgb(image::Matrix{Float64};
        color_scheme::Symbol=:none)

Version without rescaling
"""
function _to_rgb(image::Matrix{Float64};
        color_scheme::Symbol=:none)
    haskey(colorschemes,color_scheme) ? change_default_colorscheme(color_scheme) : nothing
    return get.(DefColorScheme,image)
end
function to_rgb(image::RescaledImage{Float64};color_scheme::Symbol=:none) 
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
                    color_scheme::Symbol=:none,
                    kwargs...) 
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


 """
    radial_distribution_statistics_plot(ds::DistributionStatistics;
                show_lower_bound::Bool=false,
                show_upper_bound::Bool=false,
                show_std::Bool=true,
                is_use_student::Bool=true,
                probability::Float64=0.95,
                length_scaler::Float64=1.0,
                is_centered::Bool=true,
                bound_color::Symbol=:red,
                # plot kwargs
                label::Union{AbstractString,Nothing} =nothing,
                minorgrid::Bool=true,
                gridlinewidth=2,
                title::Union{AbstractString,Nothing} ="Average temperature radial distribution",
                framestyle::Symbol = :box,
                dpi::Int=600,
                xlabel::Union{AbstractString,Nothing} = "Distance  across the sample ,mm", 
                ylabel::Union{AbstractString,Nothing}  =  "Temperature °C",
                kwargs...)

Recipe to plot the distribution statistics [`DistributionStatistics`](@ref) object.
Returns structure of `StatDataPlot` type, which has attached recipe to plot 
radial ditribution averaged value, lower and upper confidence bounds as <d> ± std
(if `show_std` is true) and confidence bounds multiplied by the Student's coefficient
(if show_lower_bound and show_upper_bound are true) calculated for the `pobability` value, 
if `length_scaler` is provied all coordinates are multiplied by this value (can be used to convert pixels to actual units)
If `is_centered` is true coordinate goes from -L/2 to +L/2 where L is the maximum of coordinates.
Other key-word arguments are the same as for the plot functions, additional keyword arguments are
transfered directly to the plot function
Usage example:

```julia
    using Plots
    plot(radial_distribution_statistics_plot(distribution_statistics))

```
This code will plot the averaged distribution statistics together with confidence bounds

"""
function radial_distribution_statistics_plot(ds::DistributionStatistics;
                show_lower_bound::Bool=false,
                show_upper_bound::Bool=false,
                show_std::Bool=true,
                is_use_student::Bool=true,
                probability::Float64=0.95,
                length_scaler::Float64=1.0,
                is_centered::Bool=true,
                bound_color::Symbol=:red,
                # plot kwargs
                label::Union{AbstractString,Nothing} =nothing,
                minorgrid::Bool=true,
                gridlinewidth=2,
                title::Union{AbstractString,Nothing} ="Average temperature radial distribution",
                framestyle::Symbol = :box,
                dpi::Int=600,
                xlabel::Union{AbstractString,Nothing} = "Distance  across the sample ,mm", 
                ylabel::Union{AbstractString,Nothing}  =  "Temperature °C",
                kwargs...)
        points_number = length(ds.coordinate)
        if is_centered || length_scaler != 1.0
            L2plot = copy(ds.coordinate)
            if is_centered
                l_center = L2plot[int_floor(points_number/2)] 
                @. L2plot -= l_center
            end
            L2plot .*=length_scaler
        else
            L2plot=ds.coordinate
        end 

        if show_lower_bound || show_upper_bound
            (lower_bound,upper_bound) = eval_bounds(ds,is_use_student=is_use_student,probability=probability)
        else
            lower_bound = nothing
            upper_bound = nothing
        end
        plot_args = merge(NamedTuple{(:label, :minorgrid,:gridlinewidth,:title,:framestyle,:dpi,:xlabel,:ylabel,:ribbon)}(
                                (label,minorgrid,gridlinewidth,title,framestyle,dpi,xlabel,ylabel,(ds.std_D,ds.std_D))
                    ),kwargs)
        return StatDataPlot(L2plot,ds.mean_D,lower_bound,upper_bound,show_std,show_lower_bound,show_upper_bound,bound_color,plot_args)
    end   
    """
    Type for custom plot recipe

    """
    struct StatDataPlot{T}
        x::T
        y::T
        lower_bound::Union{T,Nothing}
        upper_bound::Union{T,Nothing}
        show_std::Bool
        show_lower_bound::Bool
        show_upper_bound::Bool
        bound_color::Symbol
        plot_args
    end

    @recipe function f(stat_data::StatDataPlot) 
        label-->stat_data.plot_args.label
        minorgrid-->stat_data.plot_args.minorgrid
        gridlinewidth-->stat_data.plot_args.gridlinewidth
        title-->stat_data.plot_args.title
        framestyle --> stat_data.plot_args.framestyle
        dpi-->stat_data.plot_args.dpi
        xlabel-->stat_data.plot_args.xlabel
        ylabel-->stat_data.plot_args.ylabel  
        stat_data.show_std ?  ribbon := stat_data.plot_args.ribbon : nothing 
        if stat_data.show_upper_bound
            @series begin
                seriestype := :path
                primary := false
                linecolor := stat_data.bound_color
                ribbon :=nothing
                stat_data.x, stat_data.upper_bound
            end
        end
        if stat_data.show_lower_bound
             @series begin
                seriestype := :path
                primary := false
                linecolor := stat_data.bound_color
                ribbon :=nothing
                stat_data.x, stat_data.lower_bound
            end
        end

        return (stat_data.x,stat_data.y)
    end 
    """
    angular_distribution_statistics_plot(ds::DistributionStatistics;
                show_lower_bound::Bool=true,
                show_upper_bound::Bool=true,
                show_std::Bool=true,
                probability::Float64=0.95,
                is_use_student::Bool=true,
                length_scaler::Float64=1.0,
                label=nothing,
                minorgrid=true,
                gridlinewidth=2,
                title="Average temperature angular distribution",framestyle = :box,
                dpi=600,xlabel = "Angle , °", ylabel="Temperature, °C",
                kwargs...)



 The same as [`radial_distribution_statistics_plot`](@ref) but plots averaged angular distribution
"""
function angular_distribution_statistics_plot(ds::DistributionStatistics;
                show_lower_bound::Bool=true,
                show_upper_bound::Bool=true,
                probability::Float64=0.95,
                is_use_student::Bool=true,
                length_scaler::Float64=1.0,
                label=nothing,
                minorgrid=true,
                gridlinewidth=2,
                title="Average temperature angular distribution",framestyle = :box,
                dpi=600,xlabel = "Angle , °", ylabel="Temperature, °C",
                kwargs...)

                return radial_distribution_statistics_plot(ds,
                        show_lower_bound=show_lower_bound,show_upper_bound=show_upper_bound,
                        probability=probability,is_use_student=is_use_student,
                        length_scaler=length_scaler,
                        is_centered=false,
                        label=label,
                        minorgrid=minorgrid,
                        gridlinewidth=gridlinewidth,
                        title=title,framestyle = framestyle,
                        dpi=dpi,xlabel = xlabel, ylabel=ylabel,
                        kwargs...) 

    end
"""
    along_line_distribution_plot(along_line_length,along_line_distribution;
                                        length_scaler::Float64=1.0,
                                        is_centered::Bool=true,
                                        kwargs...)


Plots temperature distribution along the line `along_line_length` - coordinates,
`along_line_distribution` - values of temperature, `length_scaler` - length scaler 
(can be used to convert pixels to the actual length)
`is_centered` - the line length is converted to the coordinates with zero value in 
the centre of the `CentredObj`

"""
function along_line_distribution_plot(along_line_length,along_line_distribution;
                                        length_scaler::Float64=1.0,
                                        is_centered::Bool=true,
                                        kwargs...)

        ds = ThermovisorImages.DistributionStatistics(along_line_length,along_line_distribution)
        return radial_distribution_statistics_plot(ds,
                        show_lower_bound=false,show_upper_bound=false,
                        probability=0.0,is_use_student=false,
                        length_scaler=length_scaler,
                        is_centered=is_centered,
                        title="Temperature distribution along the line",
                        xlabel = "Temperature °C", ylabel="Distance along the  line, mm",
                        kwargs...) 
    end
