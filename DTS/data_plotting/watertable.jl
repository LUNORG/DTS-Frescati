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

query2 = """
    SELECT datetime, length, temperature
    FROM dtsold_data
    WHERE EXTRACT(YEAR FROM datetime) IN (2021, 2022)
    AND channel = 2
    AND length BETWEEN 0.0 AND 412.0;
"""
query1 = """
    SELECT datetime, length, temperature
    FROM dtsold_data
    WHERE EXTRACT(YEAR FROM datetime) IN (2021, 2022)
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

function spline_generator(query::String, conn_string::String, start_date::DateTime, end_date::DateTime, filepath::String, lim::Float64, A::Float64, process_noise::Float64, s::Int64)
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
    spline_2d = Spline2D(time_values, fiber_lengths, Kalman_matrix, kx=3, ky=3, s=s) 

    return spline_2d, fiber_lengths, time_values, reference
end

start_date = DateTime("2021-11-13")
end_date = DateTime("2022-12-14")
filepath_channel1 = "C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\pyenv\\correct_variance_dict_channel1.json"
filepath_channel2 = "C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\pyenv\\correct_variance_dict_channel2.json"
A = 1.00                         
process_noise = 1e-2 
s1_2022 = 5800
s2_2022 = 680000
s2_2024 = 670000
s1_2024 = 31000

spline_BH10, fiber_lengths2, time_values2, reference2 = spline_generator(query2, conn_string, start_date, end_date, filepath_channel2, 0.0, A, process_noise, s2_2022)
spline_BH1, fiber_lengths1, time_values1, reference1 = spline_generator(query1, conn_string, start_date, end_date, filepath_channel1, 40.0, A, process_noise, s1_2022)
l1_ = collect(range(fiber_lengths1[8], fiber_lengths1[end-2], length = length(fiber_lengths1[8:end-2])))
L1_ = l1_ .- l1_[1]
l2_ = collect(range(fiber_lengths2[30], fiber_lengths2[end-2], length = length(fiber_lengths2[30:end-2])))
L2_ = l2_ .- l2_[1]

full_time_values = collect(0:600:time_values1[end])
num_ticks = 8
tick_indices = round.(Int, LinRange(1, length(full_time_values), num_ticks))
selected_ticks = full_time_values[tick_indices]
sel_labels = [string(Date(reference1 + Second(val))) for val in selected_ticks]
ytick_indices = round.(Int, LinRange(1, length(L1_), 6))
yselected_ticks = L1_[ytick_indices]
ysel_labels = round.(reverse(yselected_ticks))

boundaries1 = find_constant_intervals(time_values1, 700.0)
boundaries2 = find_constant_intervals(time_values2, 700.0)
T_BH1 = [spline_BH1(t,l)*valid_times(t,boundaries1) for t in full_time_values, l in fiber_lengths1]
T_BH10 = [spline_BH10(t,l)*valid_times(t,boundaries2) for t in full_time_values, l in fiber_lengths2]

plotly()
pl1 = contourf(full_time_values[1:30:end], L1_, reverse(T_BH1[1:30:end, 8:end-2], dims=2)',  
    ylabel="Depth", xlabel="Time", zlims=(4.5, 16),
    legend=true, color=:thermal, yticks=(yselected_ticks, ysel_labels),
    xticks=(selected_ticks, sel_labels), size=(800,600), colorbar=false 
)
hline!([L1_[end]-43], lw=2, label=false, color=:black)
savefig("waterlevel_BH1_2024.html")

pl2 = contourf(full_time_values[1:30:end], L2_, reverse(T_BH10[1:30:end, 30:end-2], dims=2)',  
    ylabel="Depth", xlabel="Time", zlims=(-0.5, 15),
    legend=true, color=:thermal, yticks=(yselected_ticks, ysel_labels),
    xticks=(selected_ticks, sel_labels), size=(800,600), colorbar=false 
)
hline!([L2_[end]-43], lw=2, label=false, color=:black)
savefig("waterlevel_BH10_2024.html")

pl1_2022 = contourf(time_values1[1:10:end], L1_, reverse(T_BH1[1:10:end, 8:end-2], dims=2)',  
    ylabel="Depth", xlabel="Time", zlims=(4.5, 16),
    legend=true, color=:thermal, yticks=(yselected_ticks, ysel_labels),
    xticks=(selected_ticks, sel_labels), size=(800,600), colorbar=false 
)
hline!([L1_[end]-43], lw=2, label=false, color=:black)
savefig("waterlevel_BH1_2022.html")

pl2_2022 = contourf(time_values2[1:10:end], L2_, reverse(T_BH10[1:10:end, 30:end-2], dims=2)',  
    ylabel="Depth", xlabel="Time", zlims=(-0.5, 15),
    legend=true, color=:thermal, yticks=(yselected_ticks, ysel_labels),
    xticks=(selected_ticks, sel_labels), size=(800,600), colorbar=false 
)
hline!([L2_[end]-43], lw=2, label=false, color=:black)
savefig("waterlevel_BH10_2022.html")

plBH1 = plot(pl1_2022, pl1, layout = (1, 2))
plBH10 = plot(pl2_2022, pl2, layout = (1, 2))


# ---------- T shift ----------------------
# Compute the mean for each column
average_T = mean(T_BH1[:,8:end-2], dims=2)
# Create the modified matrix with T - average_T
matrix_diff = T_BH1[:,8:end-2] .- average_T

# Plot the contour
contourplot = contourf(full_time_values[1:30:end], L1_, reverse(matrix_diff[1:30:end,:],dims=2)', 
    ylabel="Depth", xlabel="Time", zlabel="Temperature", clims=(-2.5, 2.5),
    legend=true,  color=:thermal, size=(800,600), yticks=(yselected_ticks, ysel_labels),
    xticks=(selected_ticks, sel_labels) # Reduced tick positions and labels
)
savefig(contourplot, "TempShift_BH1_2022.html")

average_T = mean(T_BH10[:,30:end-2], dims=2)
# Create the modified matrix with T - average_T
matrix_diff = T_BH10[:,30:end-2] .- average_T

# Plot the contour
contourplot = contourf(full_time_values[1:30:end], L2_, reverse(matrix_diff[1:30:end,:],dims=2)', 
    ylabel="Depth", xlabel="Time", zlabel="Temperature", clims=(-2.5, 2.5),
    legend=true,  color=:thermal, size=(800,600), yticks=(yselected_ticks, ysel_labels),
    xticks=(selected_ticks, sel_labels) # Reduced tick positions and labels
)
savefig(contourplot, "TempShift_BH10_2022.html")
