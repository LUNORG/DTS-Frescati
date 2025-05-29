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
using Interpolations

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
    WHERE EXTRACT(YEAR FROM datetime) IN (2020, 2021, 2022, 20223, 2024)
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
    
    return Kalman_matrix, fiber_lengths, time_values
end

start_date = DateTime("2020-01-01")
end_date = DateTime("2024-10-14")
filepath_channel1 = "C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\pyenv\\correct_variance_dict_channel1.json"
filepath_channel2 = "C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\pyenv\\correct_variance_dict_channel2.json"
A = 1.00                         
process_noise = 1e-2 

kalman_BH10, fiber_lengths2, time_values2 = kalman_generator(query2, conn_string, start_date, end_date, filepath_channel2, 0.0, A, process_noise)
kalman_BH1, fiber_lengths1, time_values1 = kalman_generator(query1, conn_string, start_date, end_date, filepath_channel1, 40.0, A, process_noise)

column_averages2 = mean(kalman_BH10[:,30:end-2], dims=2)[:]
column_averages1 = mean(kalman_BH1[:,8:end-2], dims=2)[:] 
spline_BH_ave2 = Spline1D(time_values2, column_averages2, k=1)
spline_BH_ave1 = Spline1D(time_values1, column_averages1, k=1)
plotly()

function fft_analysis(time_values, spline_object)
    full_time_values = collect(0:600:time_values[end])
    signal = [spline_object(t) for t in full_time_values]  # Combined time series for this fiber length
    signal_fft = fft(signal)
    n = length(signal)
    dt = 600 # Time step size
    freq = fftfreq(n, 1 / dt)  # FFT frequency bins
    # Take only positive frequencies for analysis
    half_n = div(n, 2)
    freq_pos = freq[1:half_n]  # Positive frequencies
    amp = abs.(signal_fft[1:half_n])  # Magnitude
    phase = angle.(signal_fft[1:half_n])

    return signal, freq_pos, amp, phase
end
signal1, freq_pos1, amp1, phase1 = fft_analysis(time_values1, spline_BH_ave1)
signal2, freq_pos2, amp2, phase2 = fft_analysis(time_values2, spline_BH_ave2)

positions = [1/(365*24*3600); 1/(30*24*3600); 1/(15*24*3600); 1/(7*24*3600); 1/(3*24*3600); 1/(24*3600)]
labels = ["year"; "month"; "15days"; "week"; "3days"; "day"]
freq_domain1 = plot(freq_pos1, amp1, title="Spectrum", xlims=(0, 4.6e-6), ylims=(0, 30000), 
    xlabel="frequency", label="", xticks=(positions, labels))
# pl = plot(time_domain, freq_domain, layout = 2, size=(800,500))
savefig(freq_domain1, "fft_channel1.html")

freq_domain2 = plot(freq_pos2, amp2, title="Spectrum", xlims=(0, 4.6e-6), ylims=(0, 30000), 
    xlabel="frequency", label="", xticks=(positions, labels))
savefig(freq_domain2, "fft_channel2.html")

peak_indices = sortperm(amp, rev=true)[1:5]  # Get indices of top 5 peaks
peak_freqs = freq_pos[peak_indices]  # Corresponding frequencies
peak_phs = phase[peak_indices]  # Corresponding phase values
peak_amp = amp[peak_indices]

# Initialize an empty DataFrame
Fourier_df = DataFrame(PeakFrequency=Float64[], PeakPhase=Float64[], PeakAmplitude=Float64[], PeakPeriod=Float64[], PeakPeriod_day=Float64[])
peak_period = zeros(5)
peakday = zeros(5)
for i in 1:5
    peak_period[i] = peak_freqs[i] > 0 ? peak_freqs[i]^(-1) : NaN  # Avoid division by zero
    peakday[i] = peak_freqs[i] > 0 ? (peak_freqs[i]^(-1))/(60*60*24) : NaN
    push!(Fourier_df, (peak_freqs[i], peak_phs[i], peak_amp[i], peak_period[i], peakday[i]))
end

#= 
channel 2 averaged
 Row │ PeakFrequency  PeakPhase     PeakAmplitude   PeakPeriod   PeakPeriod_day
─────┼──────────────────────────────────────────────────────────────────────────
   1 │    0.0          2.42757e-17       2.45397e6  NaN                 NaN
   2 │    3.31068e-8   1.83654           3.00624e5    3.02053e7         349.599
   3 │    1.32427e-8  -2.13009           1.90525e5    7.55133e7         873.997
   4 │    1.98641e-8   1.03478           1.17987e5    5.03422e7         582.664
   5 │    4.63495e-8   2.3331       104247.0          2.15752e7         249.713

channel 1 averaged
 Row │ PeakFrequency  PeakPhase     PeakAmplitude  PeakPeriod   PeakPeriod_day
─────┼─────────────────────────────────────────────────────────────────────────
   1 │    0.0          2.66172e-17      2.29459e6  NaN                 NaN
   2 │    3.31068e-8   1.33391          1.79293e5    3.02053e7         349.599
   3 │    1.32427e-8  -2.02322          1.17551e5    7.55133e7         873.997
   4 │    1.98641e-8   0.306362     98999.8          5.03422e7         582.664
   5 │    6.62135e-9   2.0009       84021.7          1.51027e8        1747.99
=# 

# Constants
one_day_freq = 1 / 86400    # Frequency for 1 cycle per day (assuming seconds)
one_hour_freq = 1 / 3600    # Frequency for 1 cycle per hour (assuming seconds)
tolerance = 0.0000001         # Tolerance for frequency matching

# Find amplitudes near 1 day frequency
day_indices1 = findall(f -> abs(f - one_day_freq) < tolerance, freq_pos1)
day_amp1 = amp1[day_indices1]
day_indices2 = findall(f -> abs(f - one_day_freq) < tolerance, freq_pos2)
day_amp2 = amp2[day_indices2]

# Find amplitudes near 1 hour frequency
hour_indices1 = findall(f -> abs(f - one_hour_freq) < tolerance, freq_pos1)
hour_amp1 = amp1[hour_indices1]
hour_indices2 = findall(f -> abs(f - one_hour_freq) < tolerance, freq_pos2)
hour_amp2 = amp2[hour_indices2]

mean(hour_amp2) # 9.219991499089469
mean(hour_amp1) # 11.663492641391663
mean(day_amp2) # 2908.7726706886288
mean(day_amp1) # 207.26003213227045

# ----- Fourier T outdoor
complete_AirT = JSON.parsefile("C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\DTS\\stockholm_avg_temperature.json")
AirT_dict = Dict(k => v for (k, v) in complete_AirT if Date(start_date) <= Date(k) && Date(end_date) >= Date(k))
AirT_sorted_dict = string.(sort(Date.(keys(AirT_dict))))
AirT = collect(AirT_dict[i] for i in AirT_sorted_dict) # Variances stored in the dictionary

datetimes_values = DateTime.(AirT_sorted_dict)
reference = datetimes_values[1]
datetimes_seconds = Millisecond.(datetimes_values .- reference)
time_values = [datetimes_seconds[i].value /1000 for i in 1:length(datetimes_seconds)]

n = length(AirT)
dt = 86400.0 # Time step size of 1 day
freq = fftfreq(n, 1 / dt)  # Proper FFT frequency bins
fft_result = fft(AirT)  # Apply FFT
# Take only positive frequencies for analysis
half_n = div(n, 2)
freq_pos = freq[1:half_n]  # Positive frequencies
amp = abs.(fft_result[1:half_n])  # Magnit
phase = angle.(fft_result[1:half_n])

# Find dominant frequency peaks
air_peak_indices = sortperm(amp, rev=true)[1:5]  # Get indices of top 5 peaks
air_peak_freqs = freq_pos[air_peak_indices]  # Corresponding frequencies
air_peak_phs = phase[air_peak_indices]  # Corresponding phase values
air_peak_amp = amp[air_peak_indices]

positions = [1/(5*365*24*3600); 1/(365*24*3600); 1/(6*30*24*3600); 1/(3*30*24*3600)]
labels = ["5years"; "1year"; "6months"; "3months"]
freq_domain = plot(freq_pos, amp, title="Spectrum", xlims=(0, 140e-9),
    xlabel="frequency (Hz)", label="", xticks=(positions, labels))
savefig(freq_domain, "fft_outdoorT.html")

# Initialize an empty DataFrame
air_Fourier_df = DataFrame(PeakFrequency=Float64[], PeakPhase=Float64[], PeakAmplitude=Float64[], PeakPeriod=Float64[], PeakPeriod_day=Float64[])

# Populate the DataFrame with expanded values
for i in 1:5
    air_peak_period = air_peak_freqs[i] > 0 ? air_peak_freqs[i]^(-1) : NaN  # Avoid division by zero
    air_peakday = air_peak_freqs[i] > 0 ? (air_peak_freqs[i]^(-1))/(60*60*24) : NaN
    push!(air_Fourier_df, (air_peak_freqs[i], air_peak_phs[i], air_peak_amp[i], air_peak_period, air_peakday))
end

#= 
air_Fourier_df
Row  │ PeakFrequency  PeakPhase     PeakAmplitude  PeakPeriod   PeakPeriod_day
─────┼─────────────────────────────────────────────────────────────────────────
   1 │    0.0         -3.47402e-17      13946.4    NaN                  NaN    
   2 │    3.30877e-8   2.04601           8179.62     3.02227e7          349.8  
   3 │    2.64702e-8  -0.738257          2065.13     3.77784e7          437.25 
   4 │    3.97052e-8   1.89249           1339.61     2.51856e7          291.5  
   5 │    1.32351e-8  -0.608862           982.815    7.55568e7          874.5     
=#

function days_phase(ph, period)
    days_shift = ph*period/(2*pi) 
    return days_shift
end

# peak phase relative to period: 349 days
air_shift = days_phase(2.03643, 349.599) # Toutdoor 113.30779799801495
first_air = Date("2020-01-01") + Day(round(Int,giorni)) # 2020-04-23
BH1_shift = days_phase(1.33391, 349.599) # channel 1 74.2192979024725
first_BH1 = first_air + Day(round(Int,air_BH1_shift)) # 
BH10_shift = days_phase(1.83654, 349.599) # channel 2 102.18583665300274
first_BH10 = first_air + Day(round(Int,air_BH10_shift)) # 

air_BH10_shift = days_phase(2.03643-1.83654, 349.599) # 11.121961345012213
BH10_BH1_shift = days_phase(1.83654-1.33391, 349.599) # 27.966538750530223
air_BH1_shift = days_phase(2.03643-1.33391, 349.599) # 39.088500095542436

num_ticks = 8
tick_indices = round.(Int, LinRange(1, length(full_time_values), num_ticks))
selected_ticks = full_time_values[tick_indices]
sel_labels = [string(Date(start_date + Second(val))) for val in selected_ticks]

boundaries2 = find_constant_intervals(time_values2, 700.0)
boundaries1 = find_constant_intervals(time_values1, 700.0)

periodsplot = plot(full_time_values, Tair_ave, 
    xlabel = "Time", ylabel="Temperature", xticks=(selected_ticks, sel_labels), size=(800,600),
    label="External Air", legend=true)
plot!(full_time_values, [spline_BH_ave2(t)*valid_times(t, boundaries2) for t in full_time_values], label="BH10")
plot!(full_time_values, [spline_BH_ave1(t)*valid_times(t, boundaries1) for t in full_time_values], label="BH1")
vline!([air_shift*24*3600, (air_shift+349.599)*24*3600, (air_shift+2*349.599)*24*3600, (air_shift+3*349.599)*24*3600, (air_shift+4*349.599)*24*3600], lw=2, label="350-days period")
savefig(periodsplot, "periods_extAir_Temp.html")