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
    WHERE EXTRACT(YEAR FROM datetime) = 2024
    AND channel = 2
    AND length BETWEEN 0.0 AND 412.0;
"""
# For channel 1 lenghts limits to match the errors' vector size: 40.0-160.0
# For channel 2 lengths limits to match the errors' vector size: 0.0-412.0 

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

first_valid_depth = 0.0 # channel 2 -> 0.0, channel 1 -> 40.0
filepath_channel = "variance_dict_channel1.json"
filepath_channel = "variance_dict_channel2.json"
A = 1.00                         
process_noise = 1e-2 

conn = LibPQ.Connection(conn_string)
df = DataFrame(LibPQ.execute(conn, query))
LibPQ.close(conn)

df_sorted = sort(df, :datetime)
udf = unstack(df_sorted, :datetime, :length, :temperature)
udf = udf[!, 2:end]
sorted_columns = sort(parse.(Float64, names(udf)))
udf_sorted = udf[:, Symbol.(string.(sorted_columns))]
M_T = Array(udf_sorted)
temperature_matrix = coalesce.(M_T, NaN)
complete_variance_dict = JSON.parsefile(filepath)
variance_dict = Dict(k => v for (k, v) in complete_variance_dict if parse(Float64, k) > lim) # Change for channel 1(>40.0) or 2(>0.0)
variance_sorted = string.(sort(parse.(Float64, k for k in keys(variance_dict))))
R = collect(variance_dict[i] for i in variance_sorted) # Variances stored in the dictionary
n_rows = size(temperature_matrix)[1]
measurement_noise = extend_variance(R, n_rows)
Kalman_matrix = kalman_filter(temperature_matrix, process_noise, measurement_noise, A)
