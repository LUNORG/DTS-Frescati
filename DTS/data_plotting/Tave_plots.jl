using DataFrames, CSV
using LinearAlgebra
using Plots
using Dates
using LibPQ
using IniFile
using Dierckx
using LinearAlgebra
using JSON
using Statistics
using FFTW

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

query = """
    SELECT datetime, length, temperature
    FROM dtsold_data
    WHERE EXTRACT(YEAR FROM datetime) IN (2020, 2021, 2022, 2023, 2024)
    AND channel = 1
    AND length BETWEEN 40.0 AND 160.0;
"""
# For channel 1 lenghts limits to match the errors' vector size: 40.0-160.0
# For channel 2 lengths limits to match the errors' vector size: 0.0-412.0 

function find_and_remove_missing!(df::DataFrame, datetime_col::Symbol)
    missing_rows = []
    
    for row in eachrow(df)
        if any(x -> ismissing(x) || (x isa AbstractFloat && isnan(x)), row)
            push!(missing_rows, row[datetime_col])
        end
    end

    println("Datetime values with missing data: ", missing_rows)

    # Remove rows with missing values
    filter!(row -> !(row[datetime_col] in missing_rows), df)
end

function extend_variance(R_column, n_rows)
    # Replicate the column variances across all rows
    return repeat(R_column', n_rows, 1)
end

function kalman_filter(data, process_noise, measurement_noise, A)
    n_rows, n_columns = size(data)
    filtered_data = similar(data)

    # Initialize state estimates and covariances for all columns
    x_est = zeros(n_columns)
    P_est = ones(n_columns)

    for t in 1:n_rows
        for j in 1:n_columns
            # Predict Step: Includes trend adjustment
            x_pred = A * x_est[j]
            P_pred = A * P_est[j] * A + process_noise

            # Calculate Kalman Gain
            K = P_pred / (P_pred + measurement_noise[t, j])

            # Update Step
            x_est[j] = x_pred + K * (data[t, j] - x_pred)
            P_est[j] = (1 - K) * P_pred

            # Save the updated value back to the matrix
            filtered_data[t, j] = x_est[j]
        end
    end
    return filtered_data
end

function compute_error(matrix, time_values, fiber_lengths, spline2d)
    n_rows, n_cols = size(matrix)
    max_absolute_error = 0.0
    residual_total = 0.0

    for i in 1:n_rows
        for j in 1:n_cols
            # Map the matrix indices to the corresponding fiber length and time values
            y = fiber_lengths[j]
            x = time_values[i]
            
            # Evaluate the spline function at (x, y)
            predicted = spline2d(x, y)
            
            # Compute residual (signed difference)
            residual = matrix[i, j] - predicted
            residual_total += abs(residual)   # Sum of residuals to calculate mean error
           
            # Compute absolute error
            absolute_error = abs(residual)
            if absolute_error > max_absolute_error
                max_absolute_error = absolute_error  # Track maximum absolute error
            end
        end
    end
    num_elements = n_rows * n_cols 
    mean_abs_error = residual_total / num_elements  # Mean error across all elements

    return mean_abs_error, max_absolute_error
end

function Kalman_generator(query::String, conn_string::String, start_date::DateTime, end_date::DateTime, filepath::String, lim::Float64, A::Float64, process_noise::Float64)
    conn = LibPQ.Connection(conn_string)
    df = DataFrame(LibPQ.execute(conn, query))
    LibPQ.close(conn)
    df_sorted = sort(df, :datetime)
    df_sorted = filter(row -> start_date <= row.datetime <= end_date, df_sorted)
    udf = unstack(df_sorted, :datetime, :length, :temperature)
    find_and_remove_missing!(udf, :datetime)
    udf_length_sorted = udf[!, 2:end]
    sorted_columns = sort(parse.(Float64, names(udf_length_sorted)))
    udf_length_sorted_sorted = udf_length_sorted[:, Symbol.(string.(sorted_columns))]
    M_T = Array(udf_length_sorted_sorted)
    temperature_matrix = coalesce.(M_T, NaN)
    complete_variance_dict = JSON.parsefile(filepath)
    variance_dict = Dict(k => v for (k, v) in complete_variance_dict if parse(Float64, k) > lim) # Change for channel 1(>40.0) or 2(>0.0)
    variance_sorted = string.(sort(parse.(Float64, k for k in keys(variance_dict))))
    R = collect(variance_dict[i] for i in variance_sorted) # Variances stored in the dictionary
    n_rows = size(temperature_matrix)[1]
    measurement_noise = extend_variance(R, n_rows)
    Kalman_matrix = kalman_filter(temperature_matrix, process_noise, measurement_noise, A)
    fiber_lengths = coalesce(sort(unique(df.length)), NaN)
    datetimes_values = unique(udf.datetime)
    datetimes_values = collect(skipmissing(datetimes_values))
    reference = datetimes_values[1]
    datetimes_seconds = Millisecond.(datetimes_values .- reference)
    time_values = [datetimes_seconds[i].value /1000 for i in 1:length(datetimes_seconds)]
    
    return Kalman_matrix, fiber_lengths, time_values
end

start_date = DateTime("2020-01-01")
end_date = DateTime("2024-10-14")
first_valid_depth = 40.0 # channel 2 -> 0.0, channel 1 -> 40.0
filepath_channel = "C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\pyenv\\correct_variance_dict_channel1.json"
filepath_channel = "C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\pyenv\\correct_variance_dict_channel2.json"
A = 1.00                         
process_noise = 1e-2 

kalman_matrix2, fiber_lengths2, time_values2 = Kalman_generator(query, conn_string, start_date, end_date, filepath_channel, first_valid_depth, A, process_noise)
kalman_matrix1, fiber_lengths1, time_values1 = Kalman_generator(query, conn_string, start_date, end_date, filepath_channel, first_valid_depth, A, process_noise)

function find_constant_intervals(time_values::Vector{Float64}, expected_diff::Float64)
    intervals = []  # To store the (t1, t2) intervals
    n = length(time_values)

    start_index = 1  # Start of the current interval

    while start_index < n
        end_index = start_index

        # Extend the interval as long as the condition holds
        while end_index < n && time_values[end_index + 1] - time_values[end_index] < expected_diff
            end_index += 1
        end

        # If the interval has more than one valid time step, record it
        if end_index > start_index
            push!(intervals, (time_values[start_index], time_values[end_index]))
        end

        # Move to the next possible start
        start_index = end_index + 1
    end

    return intervals
end
boundaries2 = find_constant_intervals(time_values2, 700.0)
boundaries1 = find_constant_intervals(time_values1, 700.0)

function valid_times(time::Float64, intervals::Vector{Any})
    for i in 1:length(intervals)
        if (intervals[i][1]) <= time <= (intervals[i][2])
            return 1  # Return 1 if time is within any interval
        end
    end
    return NaN  # Return NaN if time does not match any interval
end

function valid_times(time::Vector{Float64}, intervals::Vector{Any})
    return map(t -> valid_times(t, intervals), time)
end

complete_AirT = JSON.parsefile("C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\DTS\\stockholm_avg_temperature.json")
AirT_dict = Dict(k => v for (k, v) in complete_AirT if Date(start_date) <= Date(k) && Date(end_date) >= Date(k))
AirT_sorted_dict = string.(sort(Date.(keys(AirT_dict))))
AirT = collect(AirT_dict[i] for i in AirT_sorted_dict) # Variances stored in the dictionary

column_averages2 = mean(kalman_matrix2, dims=2)[:]
column_averages1 = mean(kalman_matrix1, dims=2)[:] 
spline_BH_ave2 = Spline1D(time_values2, column_averages2)
spline_BH_ave1 = Spline1D(time_values1, column_averages1)

datetimes_values = DateTime.(AirT_sorted_dict)
reference = datetimes_values[1]
datetimes_seconds = Millisecond.(datetimes_values .- reference)
time_values = [datetimes_seconds[i].value /1000 for i in 1:length(datetimes_seconds)]

full_time_values = collect(0:600:time_values2[end]) # time_values2 is the shortest!

splineTair_ave = Spline1D(time_values, AirT)

Tair_ave = [splineTair_ave(x) for x in full_time_values]
BH2_ave = [spline_BH_ave2(x)*valid_times(x, boundaries2) for x in full_time_values]
BH1_ave = [spline_BH_ave1(x)*valid_times(x, boundaries1) for x in full_time_values]

num_ticks = 8
tick_indices = round.(Int, LinRange(1, length(full_time_values), num_ticks))
selected_ticks = full_time_values[tick_indices]
selected_labels = [string(Date(reference + Second(val))) for val in selected_ticks]

plotly()
pl = plot(full_time_values, Tair_ave,
    xlabel = "Time", ylabel = "Temperature",
    xticks = (selected_ticks, selected_labels),
    label = "Air aveT",
    size = (800,600))
plot!(full_time_values, BH2_ave,
    label = "Borehole10 aveT")
plot!(full_time_values, BH1_ave,
    label = "Borehole1 aveT")
peak_seconds = 118.29938377762939*24*60*60 # Using fft results
vline!([
    peak_seconds, 
    peak_seconds + 349.599*24*60*60, 
    peak_seconds + 2*349.599*24*60*60, 
    peak_seconds + 3*349.599*24*60*60, 
    peak_seconds + 4*349.599*24*60*60
], label="349days-Periods", lw=2)

savefig(pl, "periods_extAir_Temp.html")