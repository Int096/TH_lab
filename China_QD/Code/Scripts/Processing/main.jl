include("src/data_structure.jl")
include("src/helper_function.jl")
include("src/core_functions.jl")
include("src/plotter_functions.jl")

function Processing()
    # Параметры для обработки
    bin_size = 0.01 # in second
    data_folders = ["525", "525 PPD", "622 PPD"]
    ppds = ["", "with PPD", "with PPD"]
    save_folders = "$(UInt32(bin_size*1000))_ms_per_bin"

    println("####################")
    println("Size of bin: $(UInt32(bin_size*1000)) ms")
    println("Folders with raw files: $data_folders")
    println("PPDs: $ppds")
    println("Save folders: $save_folders")
    println("####################\n")
    println("Начало работы скрипта.")

    # Path 
    path_to_ChinaQD_folder = ""
    path_to_RawDatas = joinpath(
        path_to_ChinaQD_folder,
        "RawDatas"
    )
    path_to_ProcessingResults = joinpath(
        path_to_ChinaQD_folder,
        "ProcessingResults",
        save_folders
    )

    ispath(path_to_ProcessingResults) || mkdir(path_to_ProcessingResults)

    for (folder, ppd) in zip(data_folders, ppds)
        println("\nНачало работы с каталогом $folder")
        
        path_to_data_folder   = joinpath(path_to_RawDatas, folder)
        path_to_result_folder = joinpath(path_to_ProcessingResults, folder)

        ispath(path_to_result_folder) || mkdir(path_to_result_folder)

        files = readdir(path_to_data_folder)
        for file in files
            println("Обработка файла $file")

            path_to_file = joinpath(path_to_data_folder, file)

            path_to_parsed_data = joinpath(path_to_result_folder, file, "data", "parsed_data.jld2")
            path_to_binned_data = joinpath(path_to_result_folder, file, "data", "binned_data.jld2")
            path_to_FLID_data   = joinpath(path_to_result_folder, file, "data", "FLID_data.jld2")
            path_to_trace_plot  = joinpath(path_to_result_folder, file, "plots", "Trace.png")
            path_to_hist_plot   = joinpath(path_to_result_folder, file, "plots", "Hist.png")
            path_to_FLID_plot   = joinpath(path_to_result_folder, file, "plots", "FLID_h.png")

            MakeDirs(
                folder,
                path_to_ProcessingResults,
                file
            )
            ParsingData(
                path_to_file,
                path_to_parsed_data
            )
            BinningData(
                path_to_parsed_data,
                path_to_binned_data,
                bin_size
            )
            FLID_horizontal(
                path_to_parsed_data,
                path_to_FLID_data
            )

            PlotterTrace(
                path_to_binned_data,
                path_to_trace_plot,
                bin_size,
                folder,
                file,
                ppd
            )
            PlotterHist(
                path_to_binned_data,
                path_to_hist_plot,
                folder,
                file,
                ppd
            )
            PlotterFLID(
                path_to_FLID_data,
                path_to_FLID_plot,
                folder,
                file,
                ppd
            )
        end
    end
end

Processing()
