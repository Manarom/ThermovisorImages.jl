module ThermovisorImages
    using Images,ImageShow,ImageIO
    using Plots,CSV
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

    export RescaledImage,FilteredImage, MarkeredImage,
        full_image_flag,reduced_image_flag,draw!,
        draw,fit_centred_obj!,
        radius,diameter,area,side,center,
        CircleObj,SquareObj,RectangleObj,
        CentredObj,obj_from_vect,copyobj,
        fill_im!,fill_im_external!,
        filtered_mean,filtered_std,
        filter_image,filter_image!,
        marker_image,
        read_temperature_file,
        find_temperature_files,
        along_line_distribution,
        within_mask_line_distribution

    """
    ThermovisorData is a package designed to process thermal images stored as matrices.
Each  element of thermal image represents a temperature value. The package enables users to 
load images from files, calculate temperature distributions, and compute statistical analyses
for temperatures along specified lines. Thermal image can be loaded using [`read_temperature_file`](@ref).
User defined matrix can be wrapped in [`RescaledImage`](@ref) which simple maps the values to [0,1] interval. 
After that, distinct areas (relative to their surroundings,like the most heated regions within
the scene) can be fitted to regions of interest (ROIs) using [`fit_all_patterns`](@ref). 
Firther temperature distribution along the specified lined within each ROI can be obtained, and average radial and
There is also possible to evaluate the averaged statistics 
over the fitted ROI's using [`CentredObjCollectionStat`](@ref).
After that the image patterns can be labeled using the [`marker_image`](@ref) function which returns  
    
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
    const DEFAULT_FITTING_OPTIONS = Ref(Optim.Options(x_abstol=1,iterations=50))
    const DEFAULT_TEMPERATURE_FITTING_OPTIONS = Ref(Optim.Options(x_abstol=1e-4,iterations=50))

    include("CentredObj.jl") # 
    include("CentredObjStatistics.jl")
    include("ImageTypes.jl")
    include("MarkeredImage.jl")
    include("CentredObjFitters.jl")
    include("VisualizationFunctions.jl")
    

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
    function recalculate_with_new_emissivity!(image::AbstractArray,new_emissivity::Float64;
                                        image_emissivity::Float64=1.0,
                                        λ_left::Float64=14.0,
                                        λ_right::Float64=14.5,
                                        is_in_kelvins::Bool=false)
        @assert 0 < new_emissivity <= 1.0 "Wrong emissivity value emissivity must be within (0,1] region "
        
        optimizer = Optim.Brent()
        #options = DEFAULT_TEMPERATURE_FITTING_OPTIONS[]
        Threads.@sync for (i,t) in enumerate(image)
            Threads.@spawn begin
                t_meas= is_in_kelvins ? t : t + Planck.Tₖ
                t_init =  t_meas*(image_emissivity/new_emissivity)^0.25 
                func = thermal_discrepancy_fun(t_meas,λ_left,λ_right,image_emissivity,new_emissivity)
                optim_out = optimize(func,t_init/2,1.5*t_init,optimizer,rel_tol=1e-2)#,optimizer,)
                image[i] = is_in_kelvins ? optim_out.minimizer : optim_out.minimizer - Planck.Tₖ
            end
        end  
        return image
    end
    function thermal_discrepancy_fun(t_meas::T,λ_left::T,λ_right::T,e_in::T,e_out::T) where T<:Float64
        return t->(e_in*Planck.band_power(t_meas, λₗ=λ_left,  λᵣ=λ_right) - e_out*Planck.band_power(t, λₗ=λ_left,  λᵣ=λ_right))^2
    end
    """
        is_temperature_file(file_name::AbstractString)

Checks if the file with `file_name` has an appropriate name for thermovisor temperature distribution file
"""
is_temperature_file(file_name::AbstractString)=match(r"_T([1-9]|[1-9][0-9]|[1-9][0-9][0-9]|[1-9][0-9][0-9][0-9]).csv",file_name)
    
#
end
