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

query = """
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
    # mean_absolute_error, max_absolute_error = compute_error(Kalman_matrix, time_values, fiber_lengths, spline_2d)

    return spline_2d, fiber_lengths, time_values
end
start_date = DateTime("2021-11-01")
end_date = DateTime("2022-12-01")
smoothness = 10000
first_valid_depth = 40.0 # channel 2 -> 0.0, channel 1 -> 40.0
filepath_channel = "C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\pyenv\\correct_variance_dict_channel1.json"
filepath_channel = "C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\pyenv\\correct_variance_dict_channel2.json"
A = 1.00                         
process_noise = 1e-2 

spline_2d_BH2, fiber_lengths2, time_values2 = spline_generator(query, conn_string, start_date, end_date, filepath_channel, first_valid_depth, A, process_noise, smoothness)
z_values2 = [spline_2d_BH2(x, y) for x in time_values2, y in fiber_lengths2]  
dz_values2 = derivative(spline_2d_BH2, time_values2, fiber_lengths2; nux = 1, nuy=0)

T_2(t,l) = spline_2d_2(t, l)
dTdt_2(t,l) = derivative(spline_2d_2, t, l; nux = 1, nuy=0)
dTdl_2(t,l) = derivative(spline_2d_2, t, l; nux = 0, nuy=1)
l2_ = collect(range(fiber_lengths2[30], fiber_lengths2[end-2], length = length(fiber_lengths2[30:end-2])))
L2_ = l2_ .- l2_[1]

spline_2d_BH1, fiber_lengths1, time_values1 = spline_generator(query, conn_string, start_date, end_date, filepath_channel, first_valid_depth, A, process_noise, smoothness)
z_values1 = [spline_2d_BH1(x, y) for x in time_values1, y in fiber_lengths1]  
dz_values1 = derivative(spline_2d_BH1, time_values1, fiber_lengths1; nux = 1, nuy=0)

T_1(t, l) = spline_2d_1(t, l)
dTdt_1(t, l) = derivative(spline_2d_1, t, l; nux = 1, nuy=0)
dTdl_1(t,l) = derivative(spline_2d_1, t, l; nux = 0, nuy=1)
l1_ = collect(range(fiber_lengths1[7], fiber_lengths1[end-2], length = length(fiber_lengths1[7:end-2])))
L1_ = l1_ .- l1_[1]

# one frame every 1h 
gr()
t_end = 60 * 60 * 1000.0
t_ = collect(range(0, t_end, length = 1000)) # this means 1000 frames -> 1000/30s =33 sec total animation duration
anim = @animate for t in t_  # Subsample frames for better performance
    plot(T_2.(t, l1_.+18), dTdt_2.(t, l1_.+18), L1_, xlabel="Temperature", ylabel="Derivative", zlabel="Depth", size=(700,700),
    xlims=(-1.8, 8.5),
    ylims = (-8e-5, 4e-5),
    zlims = (minimum(L1_), maximum(L1_)),    
    color=:red,
    )
    plot!(T_1.(t, l1_), dTdt_1.(t, l1_), L1_,     
    color=:blue,
    )
end
gif(anim, "channel2_3d_2022.gif", fps=30)

fiber_point = 27
selected_labels = [string(Date(DateTime(start_date) + Second(val))) for val in time_values2]

animat = @animate for t in 1:100:size(z_values2, 1)
    scatter(
        [z_values2[t, fiber_point+32]], [dz_values2[t, fiber_point+32]],  # Current point
        color=:red,
        xlabel="Temperature (°C)",
        ylabel="dT/dt (°C/s)",
        xlims=(minimum(z_values2[:,fiber_point+32])-1, maximum(z_values2[:,fiber_point+32])+1),
        ylims=(minimum(dz_values2[:,fiber_point+32]), maximum(dz_values2[:,fiber_point+32])),
        title="Time Evolution length: $(fiber_lenghts2[fiber_point+32])\n Time:$(selected_labels[t])",
        legend=false,
        markersize=10,
    )
    plot!(
        z_values2[1:t, fiber_point+32], dz_values2[1:t, fiber_point+32],  # Trajectory line
        color=:red,
        linewidth=2,
    )
    scatter!(
        [z_values1[t, fiber_point+3]], [dz_values1[t, fiber_point+3]],  # Current point
        color=:blue,
        markersize=10,
    )
    plot!(
        z_values1[1:t, fiber_point+3], dz_values1[1:t, fiber_point+3],  # Trajectory line
        color=:blue,
        linewidth=2,
    )
end
gif(animat, "channels_trace_2022.gif", fps=30)

anim = @animate for t in 1:100:size(z_values2, 1)
    # Clear the previous plot
    plot()

    # Calculate the number of points to display (e.g., last 10 points)
    num_points_to_display = 600
    start_idx = max(1, t - num_points_to_display + 1)
    # println(start_idx)
    # Plot the trajectory and scatter points for z_values2
    for i in start_idx:t
        alpha_value = (i - start_idx + 1) / num_points_to_display  # Fade effect
        plot!(
            z_values2[start_idx:i, fiber_point+30], dz_values2[start_idx:i, fiber_point+30],
            color=:red, alpha=alpha_value, linewidth=2, label=""
        )
    end
    scatter!(
        [z_values2[t, fiber_point+30]], [dz_values2[t, fiber_point+30]],
        color=:red, markersize=10, label=""
    )

    # Plot the trajectory and scatter points for z_values1
    for i in start_idx:t
        alpha_value = (i - start_idx + 1) / num_points_to_display  # Fade effect
        plot!(
            z_values1[start_idx:i, fiber_point+7], dz_values1[start_idx:i, fiber_point+7],
            color=:blue, alpha=alpha_value, linewidth=2, label=""
        )
    end
    scatter!(
        [z_values1[t, fiber_point+7]], [dz_values1[t, fiber_point+7]],
        color=:blue, markersize=10, label=""
    )

    # Set plot attributes
    plot!(
        xlabel="Temperature (°C)",
        ylabel="dT/dt (°C/s)",
        xlims=(minimum(z_values2[:,fiber_point+30])-1, maximum(z_values2[:,fiber_point+30])+1),
        ylims=(minimum(dz_values2[:,fiber_point+30]), maximum(dz_values2[:,fiber_point+30])),
        title="Time Evolution length: $(fiber_lengths2[fiber_point+30]-fiber_lengths2[30])\n Time:$(selected_labels[t])",
        legend=false
    )

    savefig("frame_$desired_frame.png")
end
gif(anim, "fading_traces.gif", fps=30)

for i in 1:16
    desired_frame = div(collect(1:100:size(z_values2, 1))[end],16)*i
    t = desired_frame
    pl = plot(
        xlabel="Temperature (°C)",
        ylabel="dT/dt (°C/s)",
        xlims=(minimum(z_values2[:,fiber_point+30])-1, maximum(z_values2[:,fiber_point+30])+1),
        ylims=(minimum(dz_values2[:,fiber_point+30]), maximum(dz_values2[:,fiber_point+30])),
        title="Time Evolution length: $(L2_[fiber_point])\n Time: $(selected_labels[t])",
        legend=false
    )

    # Calculate the number of points to display
    num_points_to_display = 600
    start_idx = max(1, t - num_points_to_display + 1)

    #   Plot the trajectory and scatter points for z_values2
    for i in start_idx:t
        alpha_value = (i - start_idx + 1) / num_points_to_display  # Fade effect
        plot!(
            z_values2[start_idx:i, fiber_point+30], dz_values2[start_idx:i, fiber_point+30],
            color=:red, alpha=alpha_value, linewidth=2, label=""
        )
    end
    scatter!(
        [z_values2[t, fiber_point+30]], [dz_values2[t, fiber_point+30]],
        color=:red, markersize=10, label="" 
    )

    # Plot the trajectory and scatter points for z_values1
    for i in start_idx:t
        alpha_value = (i - start_idx + 1) / num_points_to_display  # Fade effect
        plot!(
            z_values1[start_idx:i, fiber_point+7], dz_values1[start_idx:i, fiber_point+7],
            color=:blue, alpha=alpha_value, linewidth=2, label=""
        )
    end
    scatter!(
        [z_values1[t, fiber_point+7]], [dz_values1[t, fiber_point+7]],
        color=:blue, markersize=10, label=""
    )

    # Sa ve the frame as a PNG file
    savefig(pl, "traces_sec_$desired_frame.png")
end
