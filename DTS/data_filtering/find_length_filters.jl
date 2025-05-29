using Dates
using LibPQ
using IniFile
using DataFrames
using Plots
using CSV

filename="config.ini"
ini = read(Inifile(), filename)

# Access database settings
username = get(ini, "Database", "username", "default_value")       # Access username
password = get(ini, "Database", "password", "default_value")       # Access password
host = get(ini, "Database", "host", "default_value")              # Access host
port = get(ini, "Database", "port", "default_value")   # Convert port to integer
database = get(ini, "Database", "database", "default_value")       # Access database
sslmode = get(ini, "Database", "sslmode", "default_value")         # Access SSL mode

# Create a connection string
conn_string = "host=$host port=$port dbname=$database user=$username password=$password sslmode=$sslmode"
conn = LibPQ.Connection(conn_string)

# Define the SQL query
query = """
    SELECT datetime, length, temperature
    FROM dtsold_data
    WHERE EXTRACT(YEAR FROM datetime) = 2022;
"""
# Execute the query
df = DataFrame(LibPQ.execute(conn, query))

filteredbychannel1 = filter(row -> row.channel == 1, df) 
filteredbychannel2 = filter(row -> row.channel == 2, df) 

LibPQ.close(conn)

function compute_C_vector_sum(filtered_data::DataFrame, ϵdT::Float64, Tmin::Float64, Tmax::Float64)
    df_sorted = sort(filtered_data, :datetime)
    # Unstack the sorted DataFrame (this will unstack based on sorted datetime)
    udf = unstack(df_sorted, :datetime, :length, :temperature)
    udf = udf[!, 2:end]
    sorted_columns = sort(parse.(Float64, names(udf)))
    udf_sorted = udf[:, Symbol.(string.(sorted_columns))]
    
    # Delete the last row with missing values and the first row
    #delete!(udf_length_sorted_sorted, nrow(udf_length_sorted_sorted))  
    #delete!(udf_length_sorted_sorted, 1) 

    # Convert the DataFrame to a matrix
    temperature_matrix = Matrix(udf_sorted)

    # Compute the temperature difference (T_diff)
    T_diff = (temperature_matrix[:, 3:end] .- temperature_matrix[:, 1:end-2]) ./ 2
    temperature_matrix = temperature_matrix[:, 2:end-1]

    # Condition 1: Absolute values in T_diff must be <= ϵdT
    condition1 = abs.(T_diff) .<= ϵdT
    # Condition 2: Values in temperature_matrix must be between Tmin and Tmax (inclusive)
    condition2 = (Tmin .<= temperature_matrix .<= Tmax)
    # Combine both conditions
    C = (condition1 .& condition2)
    # Sum the values of each column in C (dimension-wise sum)
    C_column_sums = sum(C, dims=1)
    # Convert the sum to a vector (flatten the result)
    C_sum_vector = vec(C_column_sums)
    
    return C_sum_vector
end

C_channel2 = compute_C_vector_sum(filteredbychannel2, 0.8, -5.0, 30.0)
C_channel1 = compute_C_vector_sum(filteredbychannel1, 0.8, -5.0, 30.0)
length_values = filteredbychannel2[:,:length]
length_values = unique(sort!(length_values))
length_values = length_values[2:end-1]

Plots.plot(length_values, C_channel2, xlabel="Length", ylabel="Sum")

# Identify contiguous regions of high values
function find_longest_plateau(indices)
    longest_start, longest_end = indices[1], indices[1]
    current_start, current_end = indices[1], indices[1]
    max_length = 0
    
    for i in 2:length(indices)
        if indices[i] == indices[i-1] + 1
            current_end = indices[i]
        else
            # Check if current plateau is the longest
            if current_end - current_start + 1 > max_length
                longest_start, longest_end = current_start, current_end
                max_length = longest_end - longest_start + 1
            end
            # Reset for next plateau
            current_start, current_end = indices[i], indices[i]
        end
    end
    
    # Final check in case the longest plateau is at the end
    if current_end - current_start + 1 > max_length
        longest_start, longest_end = current_start, current_end
    end

    return longest_start, longest_end
end

gr()

function find_length_filters(filtered_data::DataFrame, ϵdT::Float64, Tmin::Float64, Tmax::Float64)
    C_sum_vector = compute_C_vector_sum(filtered_data, ϵdT, Tmin, Tmax)
    # Threshold for identifying "high values"
    threshold = 0.9*maximum(C_sum_vector)
    # Find indices where y is above the threshold
    high_indices = findall(C_sum_vector .> threshold)

    # Find start and end indices of the longest segment of valid measurements (fiber is inside the borehole)
    start_index, end_index = find_longest_plateau(high_indices)

    length_values = df[:,:length]
    length_values = unique(sort!(length_values))
    length_values = length_values[2:end-1]

    # Retrieve the corresponding length values for the start and end
    start_length = length_values[start_index] 
    end_length = length_values[end_index]

    return start_length, end_length
end
start_length2, end_length2 = find_length_filters(filteredbychannel2, 1.0, -5.0, 30.0)
start_length1, end_length1 = find_length_filters(filteredbychannel1, 0.8, -5.0, 30.0)
# start_length: For channel2 62.957, for channel1 44.694
# end_length: For channel2 407.919, for channel1 154.271 

# Index vector normalized
C_channel2_norm = C_channel2 ./ maximum(C_channel2)

# Plotting the result
pl = plot(length_values, C_channel2_norm, label="", xlabel="Length", ylabel="Normalized Index", lw=1)
max_value = maximum(C_channel2_norm)
hline!([0.9 * max_value], label="", color=:red)
annotate!(0.8, 0.9 * max_value, text("90%", :right, 10))
vline_positions = [62.957, end_length2]  # Replace with exact positions if needed
vline!(vline_positions, linestyle=:dash, color=:red, label="")
annotate!(vline_positions[1], 0.001, text("62.957", :center, 10))
annotate!(vline_positions[2], 0.001, text("$end_length2", :center, 10))