module ThermovisorData
    using Images,ImageShow,ImageIO
    using Plots,CSV
    using Colors, ColorVectorSpace
    using Dates,Statistics,LinearAlgebra
    #using ImageSegmentation,IndirectArrays
    using Optim
    using LaTeXStrings
    using Distributions # to evaluate the Students coefficient
    using PerceptualColourMaps # to draw heatmaps
    using StaticArrays
    using Interpolations
    using  FileTypes
    using FileIO
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
    ThermovisorData

is a package designed to process thermal images stored as matrices
Each  element of thermal image represents a temperature value. The package enables users to 
load images from files, calculate temperature distributions, and compute statistical analyses
for temperatures along specified lines. It also calculates averaged angular and radial temperature
distributions (along with standard deviations) within Regions of Interest (ROIs [`CentredObj`](@ref)) 
such as  circles, squares, and rectangles. These ROI objects can be fitted to 
distinct areas (relative to their surroundings), such as the most heated regions within
the scene.
    
    """    
    ThermovisorData

    const default_images_folder = Ref(joinpath(abspath(joinpath(@__DIR__, "..")),"thermal images"))
    const FlagMatrix = Union{Matrix{Bool},BitMatrix}
    const FlagVector = Union{Vector{Bool},BitVector}
    const int_floor_abs = Int ∘ floor ∘ abs
    const int_floor = Int ∘ floor
    const int_floor_fld = Int ∘ floor ∘ fld
    const DefColorScheme = Ref("HEAT")
    const DEFAULT_FITTING_OPTIONS = Ref(Optim.Options(x_abstol=1,iterations=50))
    include("CentredObj.jl") # 
    include("CentredObjStatistics.jl")
    include("ImageTypes.jl")
    include("MarkeredImage.jl")
    include("CentredObjFitters.jl")
    include("VisualizationFunctions.jl")
    

    """
    read_temperature_file(f_name::AbstractString)

Reads temeprature file `f_name` is a full file name
"""
function read_temperature_file(f_name::AbstractString)
		if isfile(f_name)
			file_full_path = f_name
		else
			file_full_path = joinpath(default_images_folder[],f_name)
            isfile(file_full_path) ? nothing : return nothing
		end
        if Is(FileType.Image, file_full_path)
            
            im_float = Float64.(Gray.(FileIO.load(file_full_path)))
        else
            im_float = CSV.Tables.matrix(CSV.File(file_full_path,header=false,types=Float64))
        end
		pic = RescaledImage(im_float);
		creation_time = mtime(file_full_path)
		return (pic,creation_time)
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
        `is_temperature_file(file_name::AbstractString)`

Checks if the file with `file_name` has an appropriate name for thermovisor temperature distribution file
"""
is_temperature_file(file_name::AbstractString)=match(r"_T([1-9]|[1-9][0-9]|[1-9][0-9][0-9]|[1-9][0-9][0-9][0-9]).csv",file_name)
    
#
end
