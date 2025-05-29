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

function kalman_generator(query::String, conn_string::String, start_date::DateTime, end_date::DateTime, filepath::String, lim::Float64, A::Float64, process_noise::Float64)
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
    
    return Kalman_matrix, fiber_lengths, time_values, reference
end

start_date = DateTime("2020-01-01")
end_date = DateTime("2024-10-14")
filepath_channel1 = "C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\pyenv\\correct_variance_dict_channel1.json"
filepath_channel2 = "C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\pyenv\\correct_variance_dict_channel2.json"
A = 1.00                         
process_noise = 1e-2 

kalman_BH10, fiber_lengths2, time_values2, reference2 = kalman_generator(query2, conn_string, start_date, end_date, filepath_channel2, 0.0, A, process_noise)
kalman_BH1, fiber_lengths1, time_values1, reference1 = kalman_generator(query1, conn_string, start_date, end_date, filepath_channel1, 40.0, A, process_noise)

column_averages2 = mean(kalman_BH10[:,30:end-2], dims=2)[:]
column_averages1 = mean(kalman_BH1[:,7:end-2], dims=2)[:] 

complete_AirT = JSON.parsefile("C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\DTS\\stockholm_avg_temperature.json")
AirT_dict = Dict(k => v for (k, v) in complete_AirT if Date(start_date) <= Date(k) && Date(end_date) >= Date(k))
AirT_sorted_dict = string.(sort(Date.(keys(AirT_dict))))
AirT = collect(AirT_dict[i] for i in AirT_sorted_dict) # Variances stored in the dictionary

datetimes_values = DateTime.(AirT_sorted_dict)
reference = datetimes_values[1]
datetimes_seconds = Millisecond.(datetimes_values .- reference)
time_values = [datetimes_seconds[i].value /1000 for i in 1:length(datetimes_seconds)]

plotly()

using LsqFit

t2 = time_values2  # Time data
T2 = column_averages2
t1 = time_values1  # Time data
T1 = column_averages1
tair = time_values
Tair = AirT
# Define the model function
C = 2π / (365*24*3600)  # Known periodicity
model(t, p) = p[1] .+ p[2] .* sin.(C .* t .+ p[3]) 

# Initial guess for A and B
p0_2 = [mean(T2), (maximum(T2)-minimum(T2))/2, π/4]
p0_1 = [mean(T1), (maximum(T1)-minimum(T1))/2, π/4]
p0_air = [mean(Tair), (maximum(Tair)-minimum(Tair))/2, π/4]

# Perform the fitting
fit2 = curve_fit(model, t2, T2, p0_2)
fit1 = curve_fit(model, t1, T1, p0_1)
fitair = curve_fit(model, tair, Tair, p0_air)

# Extract fitted parameters
A2, B2, phi2 = fit2.param # A=7.6642645088020025, B=-4.032918896019001, phi=1.0337473191432722
A1, B1, phi1 = fit1.param # A=8.696557364776991, B=-2.24606240222075, phi=0.48833514915887644
Aair, Bair, phiAir = fitair.param # A=7.778669604641196, B=-9.768972073779782, phi=7.463181044206208
# Plot original data and fitted function
full_time_values = collect(0:600:time_values1[end])
num_ticks = 8
tick_indices = round.(Int, LinRange(1, length(full_time_values), num_ticks))
selected_ticks = full_time_values[tick_indices]
sel_labels = [string(Date(start_date + Second(val))) for val in selected_ticks]
scatter(t2, T2, label="Original Data", ylim=(1,14), markersize=2, alpha=0.5, xlabel="Time", ylabel="Temperature", xticks=(selected_ticks,sel_labels))
plot!(collect(0:600:t2[end]), model(collect(0:600:t2[end]), fit2.param), label="Fitted Sinusoid", linewidth=2, color=:red)
savefig("fitting_sin_BH10.html")
scatter(t1, T1, label="Original Data", ylim=(1,14), markersize=2, alpha=0.5, xlabel="Time", ylabel="Temperature", xticks=(selected_ticks,sel_labels))
plot!(collect(0:600:t1[end]), model(collect(0:600:t1[end]), fit1.param), label="Fitted Sinusoid", linewidth=2, color=:red)
savefig("fitting_sin_BH1.html")
scatter(tair, Tair, label="Original Data", markersize=2, alpha=0.5, xlabel="Time", ylabel="Temperature", xticks=(selected_ticks,sel_labels))
plot!(collect(0:86400:tair[end]), model(collect(0:86400:tair[end]), fitair.param), label="Fitted Sinusoid", linewidth=2, color=:red)
savefig("fitting_sin_air.html")

T2model = [model(time, fit2.param) for time in full_time_values]
T1model = [model(time, fit1.param) for time in full_time_values]
Tairmodel = [model(time, fitair.param) for time in full_time_values]

sinplot = plot(full_time_values, t->model(t, fitair.param),
    label="T outdoor", xlabel="Time", ylabel="Temperature", xticks=(selected_ticks,sel_labels))
plot!(full_time_values, t->model(t, fit2.param), label="BH10")
plot!(full_time_values, t->model(t, fit1.param), label="BH1")
savefig(sinplot, "3sin.html")

# --- Amplitude difference
abs(B2)-abs(B1) # 1.786856493798251
abs(Bair)-abs(B2) # 5.736053177760781
abs(Bair)-abs(B1) # 7.522909671559033
maximum(T2model)-minimum(T2model) # 8.06583777840654
maximum(T1model)-minimum(T1model) # 4.492124804435836
maximum(Tairmodel)-minimum(Tairmodel) # 19.537944145590398

# --- Mean value
mean(T1model) # 8.680605894191151
mean(T2model) # 7.724045625804686
mean(Tairmodel) # 7.976197450347917

# --- Phase shift ---
Δphi2_1=(phi2-phi1)
Δphi2_air=(phi2-phiAir+2π)
Δphi1_air=(phi1-phiAir+2π)
shift_BH1_BH10 = (Δphi2_1)/(C*24*3600)
shift_air_BH10 = -(Δphi2_air)/(C*24*3600)
shift_BH1_BH10 = -(Δphi1_air)/(C*24*3600)

# other method to evaluate the phase shift
tBH10_min=0.0
tBH1_min=0.0
tair_min=0.0
for (i,t) in enumerate(full_time_values)
    if T2model[i] == minimum(T2model)
        tBH10_min = t
        println("min reach after $t seconds at $i")
    end
end
# min reach after 6.1554e6 seconds
for (i,t) in enumerate(full_time_values)
    if T1model[i] == minimum(T1model)
        tBH1_min = t
        println("min reach after $t seconds at $i")
    end
end
# min reach after 8.8014e6 seconds
for (i,t) in enumerate(full_time_values[1:50343])
    if Tairmodel[i] == minimum(Tairmodel[1:50343])
        tair_min = t
        println("min reach after $t seconds")
    end
end
# min reach after 5.3136e6 seconds

shift_BH1_BH10 = (tBH1_min-tBH10_min)/(24*3600) # 31.680555555555557 days
shift_air_BH10 = (tBH10_min-tair_min)/(24*3600) # 8.5 days
shift_air_BH1 = (tBH1_min-tair_min)/(24*3600) # 40.18055555555556 days

tBH10_max=0.0
t_BH1_max=0.0
for (i,t) in enumerate(full_time_values)
    if T2model[i] == maximum(T2model)
        tBH10_max = t
        println("max reach after $t seconds at iteration $i")
    end
end
# max reach after 1.84638e7 seconds
# max reach after 4.99998e7 seconds
# max reach after 8.15358e7 seconds at full_time_values[135894]
# max reach after 1.130718e8 seconds
# max reach after 1.446078e8 seconds
for (i,t) in enumerate(full_time_values)
    if T1model[i] == maximum(T1model)
        tBH1_max = t
        println("max reach after $t seconds at iteration $i")
    end
end
# max reach after 2.13342e7 seconds
# max reach after 5.28702e7 seconds
# max reach after 8.44062e7 seconds at full_time_values[140678]
# max reach after 1.159422e8 seconds
# max reach after 1.474782e8 seconds