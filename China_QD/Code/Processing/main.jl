include("src/core.jl")
include("src/plotters.jl")

function main()
    #########################
    #  Информация о данных
    #########################
    folders_with_raw_data = ["525", "525_PPD", "622_PPD"]
    bin_size = 0.001
    save_folder = "$(UInt32(bin_size*1000))_ms_per_bin"

    println("####################")
    println("Size of bin: $(UInt32(bin_size*1000)) ms")
    println("Folders with raw files: $folders_with_raw_data")
    println("Save folders: $save_folder")
    println("####################\n")
    println("Начало работы скрипта.")

    ispath("Results") || mkdir("Results")

    for (i, folder_with_data) in enumerate(folders_with_raw_data)
        println("Обработка каталога $folder_with_data; ($i/$(length(folders_with_raw_data)))")

        raw_data_files = readdir(joinpath("RawData", folder_with_data))

        ispath("Results/"*folder_with_data) || mkdir("Results/"*folder_with_data)
        ispath("Results/"*folder_with_data*"/"*save_folder) || mkdir("Results/"*folder_with_data*"/"*save_folder)

        path_to_folder_with_results = joinpath(
            "Results", folder_with_data, save_folder
        )

        for (j, raw_data_file) in enumerate(raw_data_files)
            println("\tОбработка файла $raw_data_file; ($j/$(length(raw_data_files)))")

            # Создаем каталоги data/ и plots/ с их внутренней структурой
            ispath(path_to_folder_with_results*"/data") || mkdir(path_to_folder_with_results*"/data")
            ispath(path_to_folder_with_results*"/data/ParsedData") || mkdir(path_to_folder_with_results*"/data/ParsedData")
            ispath(path_to_folder_with_results*"/data/BinnedData") || mkdir(path_to_folder_with_results*"/data/BinnedData")
            ispath(path_to_folder_with_results*"/data/newFLID") || mkdir(path_to_folder_with_results*"/data/newFLID")

            ispath(path_to_folder_with_results*"/plots") || mkdir(path_to_folder_with_results*"/plots")
            ispath(path_to_folder_with_results*"/plots/Trace") || mkdir(path_to_folder_with_results*"/plots/Trace")
            ispath(path_to_folder_with_results*"/plots/Hist") || mkdir(path_to_folder_with_results*"/plots/Hist")
            ispath(path_to_folder_with_results*"/plots/newFLID") || mkdir(path_to_folder_with_results*"/plots/newFLID")

            path_to_file_with_raw_data = joinpath(
                "RawData", folder_with_data, raw_data_file
            )
            path_to_file_with_parsed_data = joinpath(
                path_to_folder_with_results,
                "data", "ParsedData", folder_with_data*"_"*raw_data_file*"_parsed_data.jld2"
            )
            path_to_file_with_binned_data = joinpath(
                path_to_folder_with_results,
                "data", "BinnedData", folder_with_data*"_"*raw_data_file*"_binned_data.jld2"
            )
            path_to_file_with_newFLID_data = joinpath(
                path_to_folder_with_results,
                "data", "newFLID", folder_with_data*"_"*raw_data_file*"_newFLID_data.jld2"
            )

            path_to_png_with_trace = joinpath(
                path_to_folder_with_results,
                "plots", "Trace", folder_with_data*"_"*raw_data_file*"_trace.png"
            )
            path_to_png_with_histogram = joinpath(
                path_to_folder_with_results,
                "plots", "Hist", folder_with_data*"_"*raw_data_file*"_histogram.png"
            )
            path_to_png_with_newFLID = joinpath(
                path_to_folder_with_results,
                "plots", "newFLID", folder_with_data*"_"*raw_data_file*"_newFLID.png"
            )

            println("\t\tЧтение файла.")
            ParseDataFromBinToJLD(
                path_to_file_with_raw_data,
                path_to_file_with_parsed_data
            )
            println("\t\tБинирование данных")
            BinningData(
                path_to_file_with_parsed_data,
                path_to_file_with_binned_data,
                bin_size
            )
            println("\t\tОтрисовка траектории")
            PlotterTrace(
                path_to_file_with_binned_data,
                path_to_png_with_trace,
                bin_size,
                folder_with_data,
                raw_data_file
            )
            println("\t\tОтрисовка гистограммы")
            PlotterHistogram(
                path_to_file_with_binned_data,
                path_to_png_with_histogram,
                bin_size,
                folder_with_data,
                raw_data_file
            )
            newFLIDData(
                path_to_file_with_parsed_data, 
                path_to_file_with_newFLID_data
            )
            PlotterFLID(
                path_to_file_with_newFLID_data,
                path_to_png_with_newFLID,
                folder_with_data,
                raw_data_file
            )
        end
    end
end

main()