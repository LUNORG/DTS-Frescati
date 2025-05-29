using Dates
using LibPQ
using IniFile
using DataFrames
using Plots
using Dierckx
using LinearAlgebra
using JSON
using CSV
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
conn = LibPQ.Connection(conn_string)

query = """
    SELECT timestamp, length, temperature
    FROM frescati_data
    WHERE channel = 1
"""
# For channel 1 lenghts limits to match the errors' vector size: 40.0-160.0
# For channel 2 lengths limits to match the errors' vector size: 0.0-412.0 

# Execute the query
df = DataFrame(LibPQ.execute(conn, query)) #
histogram_datetime_df = histogram(DateTime.(df.timestamp), bins=30, xlabel="Time", ylabel="Frequency", title="Histogram of Timestamps")

LibPQ.close(conn)

# Sort the DataFrame by :datetime (ascending order of rows)
df_sorted = sort(df, :timestamp)
end_date = DateTime("2024-10-14")
df_sorted = filter(row -> row.datetime <= end_date, df_sorted)
# Unstack the sorted DataFrame (this will unstack based on sorted timestamp)
udf = unstack(df, :timestamp, :length, :temperature) #UNSTACK
df_wide = unstack(df, :temperature, :length, combine=last)
T_matrix = df_wide[1:1:end, 2:10:end]
convert(Array, T_matrix)
mat_convert = Matrix(T_matrix)
heatmap(mat_convert)
# Function to find rows with missing or NaN values
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
find_and_remove_missing!(udf, :datetime)

# Extract the sorted matrix (excluding the first column which corresponds to the datetime)
udf = udf[!, 2:end]
sorted_columns = sort(parse.(Float64, names(udf)))
udf_sorted = udf[:, Symbol.(string.(sorted_columns))]
# Convert the sorted DataFrame to a matrix
M_T = Array(udf_sorted)

filepath_channel1 = "C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\pyenv\\correct_variance_dict_channel1.json"
filepath_channel2 = "C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\pyenv\\correct_variance_dict_channel2.json"
# Load the correct JSON file into a Julia dictionary
complete_variance_dict = JSON.parsefile(filepath_channel2)
variance_dict = Dict(k => v for (k, v) in complete_variance_dict if parse(Float64, k) > 0.0)
variance_sorted = string.(sort(parse.(Float64, k for k in keys(variance_dict))))

# Replace missing with NaN
temperature_matrix = coalesce.(M_T, NaN)

σ² = collect(variance_dict[i] for i in variance_sorted) # Variances stored in the dictionary
R = σ²

function extend_variance(R_column, n_rows)
    # Replicate the column variances across all rows
    return repeat(R_column', n_rows, 1)
end

n_rows, n_columns = size(temperature_matrix)
# Extend the column variances across rows
measurement_noise = extend_variance(R, n_rows)

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

# Define parameters
A = 1.00                         # Slight upward trend
process_noise = 1e-2             # Small process noise for stable data

# Apply the Kalman filter
Kalman_matrix = kalman_filter(temperature_matrix, process_noise, measurement_noise, A)

#------ Spline interpolation
fiber_lengths = sort(unique(df.length))
fiber_lenghts = coalesce.(fiber_lengths, NaN)
datetimes_values = unique(udf.datetime)
datetimes_values = collect(skipmissing(datetimes_values))
reference = datetimes_values[1]
datetimes_seconds = Millisecond.(datetimes_values .- reference)
time_values = [datetimes_seconds[i].value /1000 for i in 1:length(datetimes_seconds)]
# time_values[2:end].-time_values[1:end-1]|>unique
full_time_values = collect(0:600:time_values[end])

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
boundaries = find_constant_intervals(time_values, 700.0)

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

# Do not use k<3 otherwise the interpolation between the missing time data is very wrong
spline_2d_50000 = Spline2D(time_values, fiber_lengths, Kalman_matrix, kx=3, ky=3, s=400000) 
z_values = [spline_2d_20000(x, y)*valid_times(x, boundaries) for x in full_time_values, y in fiber_lengths]  

# trying with different values of s to estimate the best in term of smoothness and accuracy
spline_2d_25000 = Spline2D(time_values, fiber_lengths, Kalman_matrix, kx=3, ky=3, s=500000) 
spline_2d_20000 = Spline2D(time_values, fiber_lengths, Kalman_matrix, kx=3, ky=3, s=600000) 

plotly()
plot(t->spline_2d_50000(t, fiber_lenghts[50]), 0, time_values[end])
plot!(t->spline_2d_25000(t, fiber_lenghts[50]), 0, time_values[end])
plot!(t->spline_2d_20000(t, fiber_lenghts[50]), 0, time_values[end])
# raw data downsampled
scatter!(time_values, Kalman_matrix[:,50])

# now the time dimension 
plot(x->spline_2d_30000(time_values[4020], x), fiber_lenghts[1], fiber_lenghts[end])
plot!(x->spline_2d_5000(time_values[4020], x), fiber_lenghts[1], fiber_lengths[end])
plot!(x->spline_2d_2000(time_values[4020], x), fiber_lenghts[1], fiber_lenghts[end])
# raw data downsampled
scatter!(fiber_lengths, Kalman_matrix[4020,:])

mean_absolute_error, max_absolute_error = compute_error_metrics(Kalman_matrix, time_values, fiber_lengths, spline_2d_20000)

#------------- Plotting
# For channel2 we start at length = 58 --> fiber_lengths[30] and stop at 407.919m --> fiber_lenghts[202=end-2]
# For channel 1 we start at length = 54 --> fiber_lengths[8] and stop at 154.271 --> fiber_lenghts[57=end-2] 

# Number of desired ticks
num_ticks = 8
# Generate evenly spaced indices for the ticks
tick_indices = round.(Int, LinRange(1, length(full_time_values), num_ticks))
selected_ticks = full_time_values[tick_indices]
selected_labels = [string(Date(reference + Second(val))) for val in selected_ticks]
# Generate the plot with reduced xticks
# Create a 3D surface plot
surfacepl = surface(full_time_values[1:50:end], fiber_lenghts[30:end-2], z_values[1:50:end,30:end-2]', 
    ylabel="Fiber Lengths", xlabel="Time", zlabel="Temperature",
    title="Surface Plot of Temperature Spline Function", legend=false,
    xticks=(selected_ticks, selected_labels), size=(800,600)    
)
savefig(surfacepl, "surface_channel2_2024.html")

contourpl = contourf(full_time_values[1:50:end], fiber_lenghts[30:end-2], z_values[1:50:end,30:end-2]', 
    ylabel="Fiber Lengths", xlabel="Time", zlabel="Temperature",
    title="Contour Plot of Temperature Spline Function", legend=true,  color=:thermal,
    xticks=(selected_ticks, selected_labels), # Reduced tick positions and labels
    size=(800,600)
)
savefig(contourpl, "contour_channel2_2024.html")

heatmappl = heatmap(full_time_values[1:50:end], fiber_lenghts[30:end-2], z_values[1:50:end, 30:end-2]', 
    ylabel="Fiber Lengths", xlabel="Time", zlabel="Temperature",
    title="Heatmap of Temperature Spline Function", legend=true,
    xticks=(selected_ticks, selected_labels), # Reduced tick positions and labels
    size = (800,600)
)
savefig(heatmappl, "heatmap_channel2_2024.html")

#-------------ERROR EVALUATION --------------------
# Function to compute general error metrics
function compute_error_metrics(matrix, time_values, fiber_lengths, spline2d)
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

mean_absolute_error, max_absolute_error = compute_error_metrics(Kalman_matrix, time_values, fiber_lengths, spline_2d_20000)
errors = Dict(
    "mean_absolute_error" => mean_absolute_error,
    "residuals_sum" => spline_2d_20000.fp,
    "max_absolute_error" => max_absolute_error
)

# Save errors to a JSON file
open("errors_channel2_2024.json", "w") do io
    JSON.print(io, errors)
end

#------------- T-Tave 
# Compute the mean for each column
average_T = mean(z_values[:,32:end-2], dims=2)
# Create the modified matrix with T - average_T
matrix_diff = z_values[:,32:end-2] .- average_T

# Plot the contour
contourplot = contourf(full_time_values[1:50:end], fiber_lenghts[32:end-2], matrix_diff[1:50:end,:]', 
    ylabel="Fiber Lengths", xlabel="Time", zlabel="Temperature", clims=(-2.5, 2.5),
    title="Contour Plot of Temperature Spline Function", legend=true,  color=:thermal,
    xticks=(selected_ticks, selected_labels) # Reduced tick positions and labels
)
savefig(contourplot, "TempDepthShift_channel2_2024.html")

#----------- Matrix of temperature derivative over time
time_derivative_matrix = derivative(spline_2d_20000, time_values, fiber_lengths, 1, 0)
time_derivative_df = DataFrame(time_derivative_matrix, :auto)
# Save the DataFrame to a CSV file
CSV.write("channel2_2024_derivative_matrix.csv", time_derivative_df)

#-----------------MONTH PER MONTH------------------
function process_monthly_data(df_sorted::DataFrame, fiber_lengths::Vector{Union{Missing, Float64}}, process_noise::Float64, measurement_noise::Matrix{Float64}, A::Float64, smoothness::Int64)
    # Month name mapping
    month_names = Dict(1 => "jan", 2 => "feb", 3 => "mar", 4 => "apr", 5 => "may", 6 => "jun", 
                       7 => "jul", 8 => "aug", 9 => "sep", 10 => "oct", 11 => "nov", 12 => "dec")

    # Ensure datetime column is of type DateTime
    df_sorted.datetime = DateTime.(df_sorted.datetime)
    months = unique(month.(df_sorted.datetime))

    for i in 1:length(months)
        month_name = month_names[months[i]]  # Get the month name
        println("Processing data for month: ", month_name)

        # Filter the DataFrame for the current month
        monthly_df = filter(row -> month(row.datetime) == months[i], df_sorted)

        # Unstack the DataFrame to get the matrix form
        monthly_udf = unstack(monthly_df, :datetime, :length, :temperature)
        find_and_remove_missing!(monthly_udf, :datetime)

        # Extract the sorted matrix (excluding datetime column)
        monthly_data = monthly_udf[!, 2:end]
        monthly_columns = sort(parse.(Float64, names(monthly_data)))
        monthly_data_sorted = monthly_data[:, Symbol.(string.(monthly_columns))]

        # Convert the DataFrame to a matrix
        monthly_matrix = Array(monthly_data_sorted)
        monthly_matrix = coalesce.(monthly_matrix, NaN)

        # Apply Kalman filtering column-wise
        smooth_matrix = kalman_filter(monthly_matrix, process_noise, measurement_noise, A)

        # Spline interpolation
        unique_dates = unique(monthly_udf.datetime)
        unique_dates = collect(skipmissing(unique_dates))

        ref_date = unique_dates[1]
        date_diffs = Millisecond.(unique_dates .- ref_date)
        time_values = [date_diffs[i].value / 1000 for i in 1:length(date_diffs)]
        full_time_values = collect(0:600:time_values[end])

        num_ticks = 8
        # Generate evenly spaced indices for the ticks
        tick_indices = round.(Int, LinRange(1, length(full_time_values), num_ticks))
        selected_ticks = full_time_values[tick_indices]
        selected_labels = [string(Date(ref_date + Second(val))) for val in selected_ticks]        
        intervals = find_constant_intervals(time_values, 700.0)

        spline_2d = Spline2D(time_values, fiber_lengths, smooth_matrix, kx=3, ky=3, s=smoothness)
        z_values = [spline_2d(x, y)*valid_times(x, intervals) for x in full_time_values, y in fiber_lengths]  

        surfacepl = surface(full_time_values[1:5:end], fiber_lenghts[32:end-2], z_values[1:5:end, 32:end-2]', 
            ylabel="Fiber Lengths", xlabel="Time", zlabel="Temperature",
            title="Surface Plot of Temperature Spline Function", legend=false,
            xticks=(selected_ticks, selected_labels), size=(800,600) 
        )
        savefig(surfacepl, "surface_channel2_$(month_name)2024.html")
    
        contourpl = contourf(full_time_values[1:5:end], fiber_lenghts[32:end-2], z_values[1:5:end, 32:end-2]',  
            ylabel="Fiber Lengths", xlabel="Time", zlabel="Temperature",
            title="Contour Plot of Temperature Spline Function", legend=true,  color=:thermal, 
            xticks=(selected_ticks, selected_labels), size=(800,600) 
        )
        savefig(contourpl, "contour_channel2_$(month_name)2024.html")
    
        heatmappl = heatmap(full_time_values[1:5:end], fiber_lenghts[32:end-2], z_values[1:5:end, 32:end-2]',  
            ylabel="Fiber Lengths", xlabel="Time", zlabel="Temperature",
            title="Heatmap of Temperature Spline Function", legend=true,
            xticks=(selected_ticks, selected_labels), size=(800,600) 
        )
        savefig(heatmappl, "heatmap_channel2_$(month_name)2024.html")

        # Compute errors
        mean_absolute_error, max_absolute_error = compute_error_metrics(smooth_matrix, time_values, fiber_lengths, spline_2d)
        errors = Dict(
            "mean_absolute_error" => mean_absolute_error,
            "residuals_sum" => spline_2d.fp,
            "max_absolute_error" => max_absolute_error
        )

        # Save errors to JSON
        open("errors_channel2_$(month_name)2024.json", "w") do io
            JSON.print(io, errors)
        end

        # Derivative of temperature
        derivative_matrix = derivative(spline_2d, time_values, fiber_lengths, 1, 0)
        derivative_df = DataFrame(derivative_matrix, :auto)
        CSV.write("channel2_$(month_name)2024_derivative_matrix.csv", derivative_df)

        println("Month $(month_name) processing complete.")
    end
end

process_monthly_data(df_sorted, fiber_lengths, process_noise, measurement_noise, A, 20000)