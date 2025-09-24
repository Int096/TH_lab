DataReadFromFile(path::String) = hton.(reinterpret(UInt32, read(path)))

function GetCorrelation(ch1::Vector{Float64}, 
                        ch2::Vector{Float64}, 
                        nums_of_points::Int32,
                        time_between::Float64,
                        end_coef::Float64
                        )
    min_ch1, max_ch1 = minimum(ch1), maximum(ch1)
    min_ch2, max_ch2 = minimum(ch2), maximum(ch2)
    num_ch1, num_ch2 = length(ch1), length(ch2)

    experiment_time = max(max_ch1, max_ch2) - min(min_ch1, min_ch2)

    bin_list = vcat(collect(1:100:1000), round.(Int, exp10.(range(3, log10(end_coef * experiment_time / time_between), length=nums_of_points))))

    # Оставляем только нечетные
    bin_list[bin_list .> 1000 .&& iseven.(bin_list)] .+= 1
    bin_list = unique(bin_list)

    nums_of_points = length(bin_list)
    bin_times = bin_list .* time_between 

    points_between = 0.5 * (bin_times[1:end-1] + bin_times[2:end])
    left_points    = vcat(0.5 * time_between, points_between)
    right_points   = vcat(points_between, bin_times[end])

    g2 = zeros(Float64, nums_of_points)
    start_point = ones(Int, nums_of_points)
    st = fill(num_ch1, nums_of_points)
    cc = zeros(Int, nums_of_points)
    j_min = ones(Int, nums_of_points)
    j_max = ones(Int, nums_of_points)

    for ii in 1 : num_ch1
        # println(ii/num_ch1*100)
        ii_time = ch1[ii]
        for m in 1:nums_of_points
            if ii_time + left_points[m] > max_ch1
                cc[m] += 1
                if cc[m] == 1
                    st[m] = ii - 1
                end
                break
            end

            if ii == 1
                if m > 1
                    j_min[m] = j_min[m-1]
                    j_max[m] = j_max[m-1]
                else
                    j_min[m] = 1
                    j_max[m] = 1
                end 
            end

            for jj in j_min[m]:num_ch2
                if ii == 1 && ch2[jj] - min(min_ch1, min_ch2) < left_points[m]
                    start_point[m] = jj
                end
                if jj == 1
                    if ch2[jj] > ii_time + left_points[m]
                        j_min[m] = jj
                        break
                    end
                elseif ch2[jj-1] < ii_time + left_points[m] && ch2[jj] > ii_time + left_points[m]
                    j_min[m] = jj
                    break
                end
            end

            for jj in j_max[m]:num_ch2
                if jj == 1
                    if ch2[jj] > ii_time + right_points[m]
                        j_max[m] = jj
                        break
                    end
                elseif ch2[jj-1] < ii_time + right_points[m] && ch2[jj] > ii_time + right_points[m]
                    j_max[m] = jj
                    break
                end
            end

            g2[m] += j_max[m] - j_min[m]
        end
    end
    get_cor = g2 ./ st * experiment_time ./(right_points .- left_points) ./ num_ch2

    return get_cor, bin_times, g2
end