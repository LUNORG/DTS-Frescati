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
host = get(ini, "Database", "host", "default_value")               # Access host
port = get(ini, "Database", "port", "default_value")               # Convert port to integer
database = get(ini, "Database", "database", "default_value")       # Access database
sslmode = get(ini, "Database", "sslmode", "default_value")         # Access SSL mode

# Create a connection string
conn_string = "host=$host port=$port dbname=$database user=$username password=$password sslmode=$sslmode"

query2 = """
    SELECT datetime, length, temperature
    FROM dtsold_data
    WHERE EXTRACT(YEAR FROM datetime) IN (2020, 2021, 2022, 2023, 2024)
    AND channel = 2
    AND length BETWEEN 0.0 AND 412.0;
"""
query1 = """
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

function month_average(query::String, conn_string::String, start_date::DateTime, end_date::DateTime, filepath::String, lim::Float64, A::Float64, process_noise::Float64, a::Int64)
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
    datetimes_values = unique(udf.datetime)
    datetimes_values = collect(skipmissing(datetimes_values))
    column_ave = mean(Kalman_matrix[:,a:end-2], dims=2)[:]
    start_idx = 1  # Starting index for each segment
    years = [2020, 2021, 2022, 2023, 2024]
    monthly_avg = fill(NaN,60)

    for j in 1:5
        l = zeros(Int64, 12)
        # Compute the lengths of each segment
        for i in 1:12
            l[i] = count(d -> year(d) == years[j] && month(d) == i, datetimes_values)
        end
    
        # Compute means using segment ranges
        for i in 1:12
            index = (years[j]-2020)*12 + i
            if l[i] > 0  # Only compute mean if there are values in this range
                end_idx = start_idx + l[i] - 1   # Determine the last index for this segment
                monthly_avg[index] = mean(column_ave[start_idx:end_idx])
                start_idx = end_idx + 1  # Update the start index for the next segment
                # println("start index updated $start_idx, completed $(years[j]) and month $i")
            end
        end
    end

    return monthly_avg
end

function week_average(query::String, conn_string::String, start_date::DateTime, end_date::DateTime, filepath::String, lim::Float64, A::Float64, process_noise::Float64, a::Int64)
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
    datetimes_values = unique(udf.datetime)
    datetimes_values = collect(skipmissing(datetimes_values))
    column_ave = mean(Kalman_matrix[:,a:end-2], dims=2)[:]
    start_idx = 1  # Starting index for each segment
    years = [2020, 2021, 2022, 2023, 2024]
    weekly_avg = fill(NaN,240)

    for j in 1:5
        l = zeros(Int64, 48)
        week_index = 1  # Index for storing values in `l`

        for i in 1:12  # Iterate over months
            for j in 1:4  # Iterate over 4 weeks per month (assuming a fixed 4-week structure)
                # Count the number of dates that fall within the given year, month, and week range
                l[week_index] = count(d -> year(d) == years[j] && month(d) == i && day(d) > (j - 1) * 7.75 && day(d) <= j * 7.75, datetimes_values)
                week_index += 1  # Move to the next index in `l`
            end
        end

        for i in 1:48
            index = (years[j]-2020)*48 + i
            if l[i] > 0  # Only compute mean if there are values in this range
                end_idx = start_idx + l[i] - 1   # Determine the last index for this segment
                weekly_avg[index] = mean(column_ave[start_idx:end_idx])
                start_idx = end_idx + 1  # Update the start index for the next segment
                # println("start index updated $start_idx, completed $(years[j]) and month $i")
            end
        end
    end

    return weekly_avg
end

#= use this for daily averages but looping over the years
    daily_avg = fill(NaN,1827)
        
    l = zeros(Int64, 366)  # One entry per day in 2020
    day_index = 1  # Tracks position in l
        
    for month_idx in 1:12  # Iterate over months
        days_in_month = daysinmonth(Date(2024, month_idx, 1))  # Get correct number of days
            
        for day_idx in 1:days_in_month  # Iterate over days
            l[day_index] = count(d -> year(d) == years[5] && month(d) == month_idx && day(d) == day_idx, datetimes_values)
            day_index += 1  # Move to next day in l
        end
    end
                
    for i in 1:length(l)
        index = (years[5]-2020)*length(l) + i
        if l[i] > 0  # Only compute mean if there are values in this range
            end_idx = start_idx + l[i] - 1   # Determine the last index for this segment
            daily_avg[index] = mean(column_ave[start_idx:end_idx])
            start_idx = end_idx + 1  # Update the start index for the next segment
            # println("start index updated $start_idx, completed $(years[j]) and month $i")
        end
    end
=#

start_date = DateTime("2020-01-01")
end_date = DateTime("2024-10-14")
filepath_channel1 = "C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\pyenv\\correct_variance_dict_channel1.json"
filepath_channel2 = "C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\pyenv\\correct_variance_dict_channel2.json"
A = 1.00                         
process_noise = 1e-2 

month_ave_BH10 = month_average(query2, conn_string, start_date, end_date, filepath_channel2, 0.0, A, process_noise, 30)
month_ave_BH1 = month_average(query1, conn_string, start_date, end_date, filepath_channel2, 40.0, A, process_noise, 8)

# Define the months for the columns
months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

# Initialize an empty DataFrame to store the results
df_BH1 = DataFrame(Year = 2020:2024)
df_BH10 = DataFrame(Year = 2020:2024)

# Add the months columns to the DataFrame
for i in 1:12
    df_BH1[!, months[i]] = month_ave_BH1[i:12:end]
    df_BH10[!, months[i]] = month_ave_BH10[i:12:end]
end

# Write the DataFrame to a CSV file
CSV.write("monthly_Tave_BH10.csv", df_BH10)
CSV.write("monthly_Tave_BH1.csv", df_BH1)

weekly_avg_BH10 = week_average(query2, conn_string, start_date, end_date, filepath_channel2, 0.0, A, process_noise, 30)
weekly_avg_BH1 = week_average(query1, conn_string, start_date, end_date, filepath_channel2, 40.0, A, process_noise, 8)

function generate_fixed_week_labels()
    labels = [string(week) * month for month in months for week in 1:4]
    return labels
end
week_labels = generate_fixed_week_labels()

# Initialize an empty DataFrame to store the results
weekly_BH1 = DataFrame(Year = 2020:2024)
weekly_BH10 = DataFrame(Year = 2020:2024)

# Add the months columns to the DataFrame
for i in 1:48
    weekly_BH1[!, week_labels[i]] = weekly_avg_BH1[i:48:end]
    weekly_BH10[!, week_labels[i]] = weekly_avg_BH10[i:48:end]
end

# Write the DataFrame to a CSV file
CSV.write("weekly_Tave_BH10.csv", weekly_BH10)
CSV.write("weekly_Tave_BH1.csv", weekly_BH1)

# Write the DataFrame to a CSV file
CSV.write("daily_Tave_BH10.csv", DataFrame(Value = daily_avg_BH10))
CSV.write("daily_Tave_BH1.csv", DataFrame(Value = daily_avg_BH1))
