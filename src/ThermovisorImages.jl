module ThermovisorImages
    using Images,ImageShow,ImageIO
    using CSV
    #using Plots
    using RecipesBase
    using Colors, ColorVectorSpace
    using Dates,Statistics,LinearAlgebra
    #using ImageSegmentation,IndirectArrays
    using Optim
    using LaTeXStrings
    using Distributions # to evaluate the Students coefficient
    using ColorSchemes # to draw heatmaps
    using StaticArrays
    using Interpolations
    using  FileTypes # to check the type of load
    using FileIO
    using StatsBase
    import PlanckFunctions as Planck
    import ImageDraw

    """
    ThermovisorData is a package designed to process thermal images stored as matrices.
Each  element of thermal image represents a temperature value. The package enables users to 
load images from files, calculate temperature distributions, compute statistical analyses
for temperatures along specified lines and recalculate the temperature by taking into account the 
emissivity and spectral range of thermal imaging camera.  
    
    
Thermal image can be loaded using [`read_temperature_file`](@ref).
User defined matrix can be wrapped in [`RescaledImage`](@ref) which simple maps the values to [0,1] interval. 
After that, distinct areas (relative to their surroundings,like the most heated regions within
the scene) can be fitted to regions of interest (ROIs) using [`fit_all_patterns`](@ref). 
Temperature distribution along the specified lined within each ROI can be obtained, and average radial and
There is also possible to evaluate the averaged statistics 
over the fitted ROI's using [`CentredObjCollectionStat`](@ref) and recalculate the temperature for a new emissivity value
of the whole image or it's region using [`recalculate_with_new_emissivity!`](@ref)
    
    """    
    ThermovisorImages

    const default_images_folder = Ref(@eval begin 
                        images_path = joinpath(abspath(joinpath(@__DIR__(), "..")),"thermal images")
                        return isdir(images_path) ? images_path : @__DIR__()
    end)
    const FlagMatrix = Union{Matrix{Bool},BitMatrix}
    const FlagVector = Union{Vector{Bool},BitVector}
    const int_floor_abs = Int ∘ floor ∘ abs
    const int_floor = Int ∘ floor
    const int_floor_fld = Int ∘ floor ∘ fld
    const DefColorScheme = Ref(ColorSchemes.inferno)# default colorscheme 
    const DefRoiColor = Ref(RGB{Float64}(0,1,0))# default roi frame color
    const DEFAULT_FITTING_OPTIONS = Ref(Optim.Options(x_abstol=1e-1,iterations=100))
    const DEFAULT_TEMPERATURE_FITTING_OPTIONS = Ref(Optim.Options(x_abstol=1e-4,iterations=50))

    export read_temperature_file,find_temperature_files,recalculate_with_new_emissivity! 
    include("CentredObj.jl") #= export CentredObj,CircleObj,SquareObj,RectangleObj,
         c_view,line_within_mask,along_line_distribution,
         within_mask_line_points_distribution,mean_within_mask,
         along_mask_line_distribution,radial_distribution =#
    include("CentredObjStatistics.jl")#= export CentredObjCollectionStat, std_within_mask,
    radial_distribution_statistics,DistributionStatistics,
    angular_distribution_statistics,plot_radial_distribution_statistics,
    generate_random_objs=#
    include("ImageTypes.jl") #export RescaledImage,FilteredImage,filter_image,
        #filter_image!,filtered_mean,filtered_std
    include("MarkeredImage.jl") #export MarkeredImage,m_view,sort_reduce!,sort_markers!,marker_image
    include("CentredObjFitters.jl")#export fit_centred_obj!,fit_all_patterns!,fit_all_patterns
    include("VisualizationFunctions.jl")#=export draw!,draw,to_rgb,draw_line_within_mask,
    draw_line_within_mask!,change_default_colorscheme,
    change_default_roi_color=#
    

    """
    read_temperature_file(f_name::AbstractString;inverse_intensity::Bool=false)


Reads temeprature file `f_name` is a full file name, `inverse_intensity` is true if 
the intensities in the loaded file should be inverted
"""
function read_temperature_file(f_name::AbstractString;
                    inverse_intensity::Bool=false,
                    color_scheme=:none)
		if isfile(f_name)
			file_full_path = f_name
		else
			file_full_path = joinpath(default_images_folder[],f_name)
            isfile(file_full_path) ? nothing : return nothing
		end
        if Is(FileType.Image, file_full_path)
            im_data = FileIO.load(file_full_path)
            im_float = Matrix{Float64}(undef,size(im_data)[1:2])
            if haskey(colorschemes,color_scheme)
                map!(x-> getinverse(colorschemes[color_scheme], x),im_float,im_data)
            else
                if color_scheme!=:none
                    println("Entered color_scheme `$(color_scheme)` is not supported using uniform convereter")
                end
                map!( Float64∘Gray ,im_float,im_data)
            end
        else
            im_float = CSV.Tables.matrix(CSV.File(file_full_path,header=false,types=Float64))
        end
		return RescaledImage(im_float,inverse_intensity = inverse_intensity)
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
    recalculate_with_new_emissivity!(image::AbstractArray,new_emissivity::Float64,
                                        image_emissivity::Float64;
                                        λ_left::Union{Float64,Nothing}=14.0,
                                        λ_right::Union{Float64,Nothing}=14.5,
                                        is_in_kelvins::Bool=false,
                                        rel_tol::Float64=1e-3)

Function recalculates all temperatures in `image` assuming the temperatures were measured 
for the surface with emissivity `image_emissivity`; `new_emissivity` is a new value of  
emissivity, `λ_left` and `λ_right` are left and right wavelength boundaries of infrared camera spectral
range in μm. If `is_in_kelvins` is false (default) all temperatures supposed to be in Celsius and 
Kelvins otherwise.
 
To find a new temperature value `T` the function solves a non-linear equation. If both λ_left and
λ_right are not equal and none of them is `nothing`, the following equation is solved:

`new_emissivity` ⋅∫iᵦ(λ,T)dλ = `image_emissivity` ⋅∫iᵦ(λ,Tₘ)dλ

Here `iᵦ` is the Planck spectral radiance, Tₘ is the measured temperature (value of pixel intensity).
The integration is preformed over the `[λ_left,λ_right]` spectral range with the help of
`band_power` function provied by the `PlanckFunctions` package. 
    
If one of (but not both) λ_left and λ_right is set to `nothing` or both of them are equal,
than the following equaion is solved:

`new_emissivity`⋅iᵦ(λ,T) = `image_emissivity`⋅iᵦ(λ,Tₘ)
     
In both cases the equation is solved by minimizing the univariate function of least-square 
difference of left and right parts of the equation. Second version of equation is much faster,
but less precise.

"""
function recalculate_with_new_emissivity!(image::AbstractArray,new_emissivity::Float64,
                                        image_emissivity::Float64;
                                        λ_left::Union{Float64,Nothing}=14.0,
                                        λ_right::Union{Float64,Nothing}=14.5,
                                        is_in_kelvins::Bool=false,
                                        rel_tol::Float64=1e-3)
        if new_emissivity==image_emissivity
            return image 
        end                                
        @assert 0 < new_emissivity <= 1.0 && 0 < image_emissivity <= 1.0 "Wrong emissivity value emissivity must be within (0,1] region "
        @assert !(isnothing(λ_left) && isnothing(λ_right)) "At least one spectra; boundary should be settled to value"
        !isnothing(λ_left) && !isnothing(λ_right) && λ_left > λ_right && ((λ_left,λ_right) = swap(λ_left,λ_right))

        optimizer = Optim.Brent()
        #options = DEFAULT_TEMPERATURE_FITTING_OPTIONS[]
        Threads.@sync for (i,t) in enumerate(image)
            Threads.@spawn begin
                t_meas= is_in_kelvins ? t : to_K(t)
                t_init =  t_meas*(image_emissivity/new_emissivity)^0.25 
                func = thermal_discrepancy_fun(t_meas,λ_left,λ_right,image_emissivity,new_emissivity)
                optim_out = optimize(func,t_init/2,1.5*t_init,optimizer,rel_tol=rel_tol)#,optimizer,)
                image[i] = is_in_kelvins ? optim_out.minimizer : optim_out.minimizer - Planck.Tₖ
            end
        end  
        return image
    end
    function thermal_discrepancy_fun(t_meas::T,λ_left::Union{T,Nothing},λ_right::Union{T,Nothing},e_in::T,e_out::T) where T<:Float64
        
        is_single_wavelength = (isnothing(λ_left) && !isnothing(λ_right)) || (!isnothing(λ_left) && isnothing(λ_right)) || (λ_left == λ_right) 
        if is_single_wavelength
            λ_left = isnothing(λ_left) ? λ_right : λ_left
            return t->(e_in*Planck.ibb(λ_left ,t_meas) - e_out*Planck.ibb(λ_left, t))^2
        end
        return t->(e_in*Planck.band_power(t_meas, λₗ=λ_left,  λᵣ=λ_right) - e_out*Planck.band_power(t, λₗ=λ_left,  λᵣ=λ_right))^2
    end
    """
    recalculate_with_new_emissivity!(image::AbstractArray,c::CentredObj,new_emissivity::Float64;
                                                kwargs...)

Recalculates temperarture of each pixel within the `CentredObj` with `new_emissivity`
see [`recalculate_with_new_emissivity!`](@ref) for keyword arguments
"""
function recalculate_with_new_emissivity!(image::AbstractArray,c::CentredObj,new_emissivity::Float64, image_emissivity::Float64;
                                                kwargs...)
        im_view = c_view(image,c)
        recalculate_with_new_emissivity!(im_view,new_emissivity,image_emissivity;kwargs...)
        return image
    end
    """
    recalculate_with_new_emissivity!(image::AbstractArray,marker::MarkeredImage,
        label::Int,new_emissivity::Float64,image_emissivity::Float64;
        kwargs...)

Recalculates the temperarture of each pixel within image pattern `label`
of `MarkeredImage` image `marker` with `new_emissivity` assuming `image_emissivity`
be the emissvity settled during measurements. Both `marker` and
`image` should be of the same size. See [`recalculate_with_new_emissivity!`](@ref)

"""
function recalculate_with_new_emissivity!(image::AbstractArray,marker::MarkeredImage,
        label::Int,new_emissivity::Float64,image_emissivity::Float64;
                                                        kwargs...)
        @assert size(image) == size(marker) "Both `image` and `marker` should be of the same size"
        im_view = m_view(image,marker,label)
        recalculate_with_new_emissivity!(im_view,new_emissivity,image_emissivity;kwargs...)
        return image
    end    
    """
        is_temperature_file(file_name::AbstractString)

Checks if the file with `file_name` has an appropriate name for thermovisor temperature distribution file
"""
is_temperature_file(file_name::AbstractString)=match(r"_T([1-9]|[1-9][0-9]|[1-9][0-9][0-9]|[1-9][0-9][0-9][0-9]).csv",file_name)
    
"""
    convert_temperature_to_emissivity(image::AbstractMatrix,
                                            real_temperature::Float64,
                                            image_emissivity::Float64;
                                            λ_left::Union{Float64,Nothing}=14.0,
                            λ_right::Union{Float64,Nothing}=14.5,
                            is_in_kelvins::Bool=false)

Evaluates the emissivity, assuming the whole `image` to be of the same temperature
`real_temperature`.
If both λ_left and λ_right are not equal and none of them is `nothing` the emissivity `ϵ` is calculated using:

ϵ =  `image_emissivity`⋅∫iᵦ(λ,Tᵣ)dλ/∫iᵦ(λ,Tₘ)dλ

Here `iᵦ` is the Planck spectral radiance, Tₘ is the measured temperature (value of pixel temperature),
Tᵣ is the real temperature of the surface.
The integration is preformed over the `[λ_left,λ_right]` spectral range with the help of
`band_power` function provied by the `PlanckFunctions` package. 
    
If one of (but not both) λ_left and λ_right is set to `nothing` or both of them are equal,
than the following equaion is solved:

ϵ =  `image_emissivity`⋅iᵦ(λ,Tᵣ)/iᵦ(λ,Tₘ)

"""
function convert_temperature_to_emissivity(image::AbstractMatrix,
                                            real_temperature::Float64,
                                            image_emissivity::Float64;
                                            λ_left::Union{Float64,Nothing}=14.0,
                            λ_right::Union{Float64,Nothing}=14.5,
                            is_in_kelvins::Bool=false)
         e_out = Matrix{Float64}(undef,size(image)...)
         fun_to = to_emissivity_fun(real_temperature,image_emissivity,λ_left,λ_right)
         Threads.@threads for (ii,t_meas) in collect(enumerate(image))
            Threads.@spawn begin
                t_meas = is_in_kelvins ? t_meas : to_K(t_meas)
                e_out[ii] = fun_to(t_meas)
            end
         end
         return e_out
         
    end
    to_K(t) = t + Planck.Tₖ
    function to_emissivity_fun(t_meas::Float64,
                            e_in::Float64,
                            λ_left::Union{Float64,Nothing},
                            λ_right::Union{Float64,Nothing})

        is_single_wavelength = (isnothing(λ_left) && !isnothing(λ_right)) || (!isnothing(λ_left) && isnothing(λ_right)) || (λ_left == λ_right) 
        if is_single_wavelength
            λ_left = isnothing(λ_left) ? λ_right : λ_left
            return t->e_in*Planck.ibb(λ_left ,t_meas)/Planck.ibb(λ_left, t)
        end
        return t->e_in*Planck.band_power(t_meas, λₗ=λ_left,  λᵣ=λ_right)/Planck.band_power(t, λₗ=λ_left,  λᵣ=λ_right)
    end
end
