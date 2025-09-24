using JLD2

begin
    parsed_data = JLD2.load_object(
        "ProcessingResults/1_ms_per_bin/525/66/data/parsed_data.jld2"
    )

    intensity = parsed_data.intensity

    
end