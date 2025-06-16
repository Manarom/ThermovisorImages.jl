using ThermovisorImages
using Test
test_data_folder = joinpath(".","test","tests data")
@testset "ThermovisorImages.jl" begin
    # csv file reading 
    rescaled_csv = read_temperature_file(joinpath(test_data_folder,"T500.csv"))
    markered_csv = marker_image(rescaled_csv)
    circle = [CircleObj()]
    fit_all_patterns!(circle,markered_csv)
    filter_image(rescaled_csv,circle[])
    @test mean_within_mask(rescaled_csv.initial,circle[]) ≈ 450 atol=20.0
    square = [SquareObj()]
    fit_all_patterns!(square,markered_csv)
    rect = [RectangleObj()]
    fit_all_patterns!(rect,markered_csv)
    @test ThermovisorImages.side(square[]) ≈ ThermovisorImages.side(rect[]) rtol=5e-2
    filtered_image_csv = filter_image(rescaled_csv,label=1)
    fit_centred_obj!(circle[],filtered_image_csv)
    fit_centred_obj!(square[],filtered_image_csv)

    @test ThermovisorImages.side(circle[]) ≈ ThermovisorImages.side(square[]) atol=35
    #image reading
    rescaled_jpg = read_temperature_file(joinpath(test_data_folder,"coins.jpg"), inverse_intensity=true)  
    markered_coins = marker_image(rescaled_jpg) 
    @test length(markered_coins)==24
    rois_coins = fit_all_patterns(rescaled_jpg)
    @test length(markered_coins)==length(rois_coins)

end
#=
test_data_folder = abspath(test_data_folder)
rescaled_csv = read_temperature_file(joinpath(test_data_folder,"T500.csv"))
rescaled_jpg = read_temperature_file(joinpath(test_data_folder,"coins.jpg"), inverse_intensity=true)  
markered_coins = marker_image(rescaled_jpg) 
rois_coins = fit_all_patterns(rescaled_jpg)


circle = [CircleObj()]
markered_csv = marker_image(rescaled_csv)
fit_all_patterns!(circle,markered_csv)
filter_image(rescaled_csv,circle[])
mean_within_mask(rescaled_csv.initial,circle[])
square = [SquareObj()]
fit_all_patterns!(square,markered_csv)
rgb_csv = ThermovisorImages.draw(rescaled_csv)
ThermovisorImages.draw!(rgb_csv,circle[])
ThermovisorImages.draw!(rgb_csv,square[])
square2 = fit_all_patterns(rescaled_csv,SquareObj)
ThermovisorImages.draw!(rgb_csv,square2[])
ThermovisorImages.draw(Float64.(markered_csv.markers))
=#