#using CairoMakie
using BenchmarkTools
using Statistics

include("QD_processing_files/QD_data_sturcture.jl")
include("QD_processing_files/QD_helper_function.jl")
include("QD_processing_files/QD_plotter_function.jl")
include("QD_processing_files/QD_processing_function.jl")

function DataProcessing(folder::String, file::String)
    path = joinpath("source_data", folder, file)
    save_path = joinpath("results")

    println("----------------------------")
    println("Начало обработки файла $path")
    println("----------------------------")

    bin_size = 0.001

    source_data = DataReadFromFile(path)
    println("Данные прочитаны")
    parsed_data = ParsingData(source_data)
    println("Данные спарсены")
    binned_data = BinningData(parsed_data.arrivals, bin_size)
    println("Данные разбинированы")
    #corr_data = CalculateCorrelateFunction(parsed_data, binned_data.intensity)
    println("Корреляционная функция вычислена")

    println("----------------------------")
    println("Конец обработки файла $path")
    println("----------------------------")
    PlotterTrace(save_path, binned_data.intensity, bin_size)
    #PlotterCorrFunction(corr_data)

end

DataProcessing("tests", "63")