using JLD2, FHist

struct WeightedTau
    taus::Float64
    weight::UInt32
end

function GenData(folders)
	for (i, folder) in enumerate(folders)
		println("($i/$(length(folders))): folder = $folder")

		path = joinpath(
			"/", "home", "arzamas", "Laboratory_TeorChem", 
			"China_QD", "Results", folder, "1_ms_per_bin", "data")
		binned_data_f = readdir(path*"/BinnedData")
		parsed_data_f = readdir(path*"/ParsedData")

		taus_dict = Dict{String, Vector{WeightedTau}}()

		for (j, (bd_file, pd_file)) in enumerate(zip(binned_data_f, parsed_data_f))
			number_of_dot_bd = split(bd_file, "_")[2]
			number_of_dot_pd = split(pd_file, "_")[2]
			if number_of_dot_bd != number_of_dot_pd
				println("ВСЕ ПРОПАЛО!!! Файлы не состыковались")
				return 
			end

			println("\t($j/$(length(binned_data_f))): $number_of_dot_bd-$number_of_dot_pd file")
			
			parsed_data = JLD2.load_object(path*"/ParsedData/"*pd_file)
			binned_data = JLD2.load_object(path*"/BinnedData/"*bd_file)

			delays = parsed_data.delays * 1e9

			intensity = binned_data.intensity
			index_max = binned_data.index_max
			index_min = binned_data.index_min

			num_of_levels = 100
			num_of_photons_per_fit = 1000

			intensity_border = range(
				start  = minimum(intensity),
				stop   = maximum(intensity),
				length = num_of_levels + 1
			)

			#println(intensity_border)
			weighted_taus = Vector{WeightedTau}()

			for m in 1:num_of_levels
				# Непосредственно разбиение по уровням
				if m == 1
					intensity_in_borders = findall(
						intensity .<= intensity_border[m+1]
					)
				else
					intensity_in_borders = findall(
						(intensity .<= intensity_border[m+1])
						.&
						(intensity .> intensity_border[m])
					)
				end

				# Поиск всех индексов
				indexes_of_photons = []
				for idx in intensity_in_borders 
					indexes_of_photons = vcat(
						indexes_of_photons, 
						index_min[idx] : index_max[idx]
					)
				end

				num_of_fit = div(length(indexes_of_photons), num_of_photons_per_fit)

				for k in 1:num_of_fit
					if k != num_of_fit
						data_idx = indexes_of_photons[
							((k-1)*num_of_photons_per_fit+1):k*num_of_photons_per_fit
						]
						data = delays[data_idx]
					else
						if k == 1
							data_idx = indexes_of_photons[
								1:end
							]
						else
							data_idx = indexes_of_photons[
								((k-1)*num_of_photons_per_fit):end
							]	
						end
						data = delays[data_idx]
					end

					hist_obj = Hist1D(
						data;
						binedges = 0.5 : (ceil(maximum(data))-0.5)
					) |> normalize
				end
			end
		end
	end
end

GenData(["800"])
