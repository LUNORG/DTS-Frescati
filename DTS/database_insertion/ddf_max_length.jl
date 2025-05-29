using DataFrames

#IF: Filtering data based on the input of our choice.
# Starting from 150m (a position reasonably part of the fiber in the borehole) the max length limit was found considering absolute temperature limits -15°C<T<40°C
function filter_dataframe(df::DataFrame, min_len::Int = 150, temp_high::Float64 = 40.0, temp_low::Float64 = -15.0)
    # Find the first row where length is at least `min_len`
    start_row = findfirst(x -> x >= min_len, df[!, :length])


    
    # If no row meets the minimum length condition, return nothing 
    if start_row === nothing
        println("No rows found with length at least $min_len.")
        return nothing
    end
    
    # Initialize the stop row to the start row, then move backwards
    #IF: searching for temperatures outside the input parametres we chose from above function, by going through each row
    stop_row = start_row
    for row in start_row:nrow(df)  #IF: df = dataframe, innehåller tabell med temperatur och längd
        temp = df[row, :temperature]
        if temp > temp_high || temp < temp_low # IF check if temeperatur is outside the chosen parametres
            stop_row = row 
            break
        end
    end
    
    # Return the sliced DataFrame from stop_row to start_row
    return df[start_row:stop_row, :]
end

# Function to find the length where temperature jump occurs
function find_length_for_temperature_jump(filtered_df::DataFrame)
    # Loop through rows starting from the second one to check temperature differences
    for i in 2:nrow(filtered_df)
        temp_diff = abs(filtered_df[i, :temperature] - filtered_df[i - 1, :temperature])
        if temp_diff >= 3  # If the temperature difference is 3 or more
            return filtered_df[i, :length]  # Return the length that meets the criteria
        end
    end
    return nothing  # Return nothing if no such length is found
end

# Function to parse a single .ddf file
function parse_ddf_file(filepath, max_length_channel_1::Ref{Union{Nothing, Float64}}, max_length_channel_2::Ref{Union{Nothing, Float64}})
    try
        lines = readlines(filepath)
        split_lines = split(lines[1], r"[\t\r]")  # Split first line to extract metadata
        metadata_lines = []
        
        for i in 1:length(split_lines)
            if split_lines[i] == "length (m)"
                break  # Stop when "length (m)" is found
            end
            push!(metadata_lines, split_lines[i])
        end
        
        # Create metadata DataFrame
        metadata_columns = metadata_lines[1:2:end]
        values = metadata_lines[2:2:end]
        metadata_df = DataFrame(Dict(metadata_columns .=> values))
        
        # Process data lines
        data = []
        # Determine the channel from metadata
        channel_forward_value = metadata_df[!, "forward channel"][1]
        channel = channel_forward_value == "channel 1" ? 1 :
                channel_forward_value == "channel 2" ? 2 : NaN
        
        for line in lines[2:end]  # Skip the first line with metadata
            values = split(line)
            if length(values) == 4  # Ensure there are exactly 4 columns
                try
                    length_val = parse(Float64, values[1])
                    temperature = parse(Float64, values[2])
                    push!(data, (length_val, temperature))  # Store data
                catch parse_error
                    println("Error parsing values in data line: $line")
                    continue  # Skip lines that cannot be parsed
                end
            end
        end
        
        # Create data DataFrame
        data_df = DataFrame(data, [:length, :temperature])
        filtered_data_df = filter_dataframe(data_df)
        # Find length that meets the condition and update max_length if needed
        length_value = find_length_for_temperature_jump(filtered_data_df)
        if channel == 1 && length_value != nothing && (isnothing(max_length_channel_1[]) || length_value > max_length_channel_1[])
            max_length_channel_1[] = length_value  # Update the max_length if a new one is found
        elseif channel == 2 && length_value != nothing && (isnothing(max_length_channel_2[]) || length_value > max_length_channel_2[])
            max_length_channel_2[] = length_value
        end
    catch e
        println("Error parsing DTS file: $e")
        return nothing  # Return nothing if any error occurs
    end
end

# Function to crawl through a folder and parse all .ddf files
function parse_all_ddf_in_folder(folder_path)
    max_length_channel_1 = Ref{Union{Nothing, Float64}}(nothing)  # Initialize max_length as nothing
    max_length_channel_2 = Ref{Union{Nothing, Float64}}(nothing)
    for (root, dirs, files) in walkdir(folder_path)
        for file in files
            if endswith(file, ".ddf")  # Check if the file is a .ddf file
                file_path = joinpath(root, file)
                parse_ddf_file(file_path, max_length_channel_1, max_length_channel_2)
            end
        end
    end
    return max_length_channel_1[], max_length_channel_2[]  # Return the max length found across all files
end

max_length = parse_all_ddf_in_folder("DTS_data\\full data set\\KTH_LiL_BH10_12")
println(max_length)
# max_length_channel_1 = 1059.289
# max_length_channel_2 = 511.408

# if i just consider KTH_LiL_BH10_12 is (158.329, 411.978) but I kept the larger limits