
include("ThermovisorImages.jl")
using .ThermovisorImages
test_data_folder = joinpath(".","test","tests data")

rescaled_csv = read_temperature_file(joinpath(test_data_folder,"T500.csv"))
markered_csv = marker_image(rescaled_csv)
circle = [CircleObj()]
fit_all_patterns!(circle,markered_csv)
new_line =200# minimum(side(fitted_obj))-5
ang_range = 0.0:1:180 # range of angles 
(R,D) = radial_distribution(rescaled_csv.initial,circle[],ang_range,line_length=new_line)

DS = ThermovisorImages.radial_distribution_statistics(R,D)

p_radial = ThermovisorImages.plot_radial_distribution_statistics(DS,show_lower_bound=true,show_upper_bound=true,probability=0.99)


p_radial = ThermovisorImages.recipe_plot_radial_distribution_statistics(DS,show_lower_bound=true,show_upper_bound=true,probability=0.99)
using Plots
plot(p_radial)