using DuckDB
using DataFrames
using Parquet
using Plots
using Dates
using Statistics
using StatsPlots
using JSON
using Measures

parquet_filename = "C:\\Users\\zezro\\Desktop\\MJ146X\\local_db.parquet"

# Execute the query and convert the result into a DataFrame
#LN: df = datafile created based on filters specified in query
# import Pkg
# Pkg.instantiate()
# this can be packed into a function

#Kalman Filter parameters
#-------------------------------------------------------------------------------------------------
A = 1.00                         # Slight upward trend (no?)
process_noise = 1e-2             # Small process noise for stable data
filepath_channel2 = "C:\\Users\\zezro\\Documents\\GitHub\\lil-test-rig-control\\DTS\\uncertanty_evaluation\\variance_dict_channel2.json" # Load the correct JSON file into a Julia dictionary

complete_variance_dict = JSON.parsefile(filepath_channel2)
variance_dict = Dict(k => v for (k, v) in complete_variance_dict if parse(Float64, k) > 0.0)
variance_sorted = string.(sort(parse.(Float64, k for k in keys(variance_dict))))
σ² = collect(variance_dict[i] for i in variance_sorted) # Variances stored in the dictionary
R = σ²
#-------------------------------------------------------------------------------------------------

#BH visaualisation analysis timeframe
#-------------------------------------------------------------------------------------------------
#min_date_analysis = DateTime("2014-01-01")
min_date_analysis = DateTime("2014-12-01")  
#max_date_analysis = DateTime("2020-12-31")
max_date_analysis = DateTime("2020-01-01")
#-------------------------------------------------------------------------------------------------

function get_temps_long(;
    channel::Int=1,
    min_length = 1520,
    max_length = 1971,
    min_timestamp::DateTime = DateTime("2014-12-01"),
    max_timestamp::DateTime = DateTime("2020-01-01"),
    parquet_filename = parquet_filename
    ) #LN: previously: parquet_filename ="local_db.parquet"

    conn = DuckDB.DB()

    min_timestamp = datetime2unix(min_timestamp)
    max_timestamp = datetime2unix(max_timestamp)

        query_filtered = """
            SELECT channel, length, timestamp, temperature
            FROM parquet_scan('$parquet_filename')
            WHERE channel = $channel
            AND length BETWEEN $min_length AND $max_length 
            AND timestamp BETWEEN $min_timestamp AND $max_timestamp
            ORDER BY timestamp;
    """

    df_filtered = DuckDB.query(conn, query_filtered) |> DataFrame
    close(conn)
    df_filtered.timestamp .= unix2datetime.(df_filtered.timestamp)
    
    return df_filtered
end
#df_l = @time get_temps_long() #LN: df_l = datafile in long format

function kalman_filter(data, process_noise, measurement_noise, A) #IF: A = avvikelse skalär
    n_rows, n_columns = size(data) #IF:rows and columns in data. 
    filtered_data = similar(data) #IF: matex to store similar results

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

function get_temps_wide(; #LN: Also gets vector of Lengths, vector of time
    channel::Int=1,
    min_length = 0.0,
    max_length = 100.0,
    min_timestamp::DateTime = DateTime("2014-12-01"),
    max_timestamp::DateTime = DateTime("2020-01-01"),
    parquet_filename = parquet_filename
    #LN: previously: parquet_filename ="local_db.parquet"

    #pass the arguments to the other function as is
    )
    df_l = get_temps_long(;
        channel = channel,
        min_length = min_length,
        max_length = max_length,
        min_timestamp = min_timestamp,
        max_timestamp = max_timestamp,
        parquet_filename = parquet_filename
    )

    df_w = unstack(df_l, :length, :temperature, combine=last) #LN: df_w = datafile in wide format

    mT = df_w[:, 3:end] |> Matrix                     #LN: mT = matrix of Temperatures
    vt = collect(df_w.timestamp)                      #LN: vt = vector of times
    vL = collect(parse.(Float64, names(df_w)[3:end])) #LN: vL = vector of Lengths
    
    return mT, vt, vL
end

function avg_rows(matrix)
    return mean(matrix, dims=1)
end

function scatterplot(channel)
    # looking at only one row, trying to identify the limits of the BHs
    mT_id_1, vt_id_1, vL_id_1 = get_temps_wide(
        channel = channel,
        min_length = -1000000,
        max_length = 1000000,
        min_timestamp = DateTime("2017-04-01"),
        max_timestamp = DateTime("2017-05-01")
    )
    #size(mT_id_1)
    #vt_id_1
    avg_mT = vec(avg_rows(mT_id_1))
     #LN: Interactable plot; zooming etc.
    gr() #LN: High res plot
    scatter(vL_id_1, avg_mT, xlabel = "Length (m)", ylabel = "Temperature (°C)", title = "Channel $channel", size=(1240, 900), margin=10mm, legend = false)
    #more advanced with go and return
    #limits_bh = Dict(12 => (1, [(023.6, 145.7); (623.6, 745.7)]))
    #simpler version
end

function normalplot(channel)
    # looking at only one row, trying to identify the limits of the BHs
    mT_id_1, vt_id_1, vL_id_1 = get_temps_wide(
        channel = channel,
        min_length = -1000000,
        max_length = 1000000,
        min_timestamp = DateTime("2017-04-01"),
        max_timestamp = DateTime("2017-05-01")
    )
    #size(mT_id_1)
    #vt_id_1
    avg_mT = vec(avg_rows(mT_id_1))
     #LN: Interactable plot; zooming etc.
    gr() #LN: High res plot
    plot(vL_id_1, avg_mT, xlabel = "Length (m)", ylabel = "Temperature (°C)", title = "Channel $channel", size=(1240, 900), margin=10mm, legend = false)
    #more advanced with go and return
    #limits_bh = Dict(12 => (1, [(023.6, 145.7); (623.6, 745.7)]))
    #simpler version
end

function extend_variance(R_column, n_rows)
    # Replicate the column variances across all rows
    return repeat(R_column', n_rows, 1)
end

#LN: Makes a mirrored depth vector, based on the length vector, for plotting.
function mirror_bh(length_vector, start_length, end_length)
    y_vec = []
    for i in length_vector
        if i-start_length<(end_length-start_length)/2
            y=start_length-i
        else
            y=i-end_length
        end
        push!(y_vec, y)
    end
    return y_vec
end

function format(dt)
    return 2*dt
end

function datewithoutyear(dt)
    res = dt - Year(dt)
    res
end

function superpos_plot(;
    channel::Int=1,
    min_length = 1520,
    max_length = 1971,
    min_timestamp::DateTime = DateTime("2014-01-01"),
    max_timestamp::DateTime = DateTime("2020-06-01"),
    parquet_filename = parquet_filename
    )

    df = get_temps_long(; 
    channel,
    min_length,
    max_length,
    min_timestamp,
    max_timestamp,
    parquet_filename
    )

    groups_time = groupby(df, :timestamp)
    df_mean = combine(groups_time, :temperature => mean)

    format(dt) = Dates.format(dt,"mm:dd") # formatting helper function

    superpos_plot = @df df_mean plot(datewithoutyear.(:timestamp), :temperature_mean, group = year.(:timestamp),ylabel="Temperature  (°C)", xlabel = "time over the year", size=(1240, 900), margin=10mm)
    vline!([DateTime("0000-03-01")], color=:gray, linestyle=:dash, label="Winter ends")
    vline!([DateTime("0000-06-01")], color=:black, linestyle=:dash, label="Spring ends")
    vline!([DateTime("0000-09-01")], color=:blue, linestyle=:dash, label="Summer ends")
    vline!([DateTime("0000-12-01")], color=:red, linestyle=:dash, label="Autumn ends")
    display(superpos_plot)

    return superpos_plot
    #savefig(superpos_plot, "C:\\Users\\zezro\\Desktop\\MJ146X\\Plots\\Superpositioned_Plot.png")
end

#LN: Makes scatterplots and normalplots with averaged temperatures for each of the four channels (gr()).
   #ch1 = scatterplot(1); savefig(ch1, "C:\\Users\\zezro\\Desktop\\MJ146X\\Plots\\Scatter_Tavg_Ch1.png")
   #ch2 = scatterplot(2); savefig(ch2, "C:\\Users\\zezro\\Desktop\\MJ146X\\Plots\\Scatter_Tavg_Ch2.png")
   #ch3 = scatterplot(3); savefig(ch3, "C:\\Users\\zezro\\Desktop\\MJ146X\\Plots\\Scatter_Tavg_Ch3.png")
   #ch4 = scatterplot(4); savefig(ch4, "C:\\Users\\zezro\\Desktop\\MJ146X\\Plots\\Scatter_Tavg_Ch4.png")

   #ch1 = normalplot(1); savefig(ch1, "C:\\Users\\zezro\\Desktop\\MJ146X\\Plots\\Plot_Tavg_Ch1.png")
   #ch2 = normalplot(2); savefig(ch2, "C:\\Users\\zezro\\Desktop\\MJ146X\\Plots\\Plot_Tavg_Ch2.png")
   #ch3 = normalplot(3); savefig(ch3, "C:\\Users\\zezro\\Desktop\\MJ146X\\Plots\\Plot_Tavg_Ch3.png")
   #ch4 = normalplot(4); savefig(ch4, "C:\\Users\\zezro\\Desktop\\MJ146X\\Plots\\Plot_Tavg_Ch4.png")



limits_bh = Dict( #LN: Assume patch 0 for Ch1
                  #LN: Comments after each are the back part of the DTS going through the BH up and down a 2nd time.
    105 => (1, (174, 626)), #105b is all kinds of bananas.
    45 => (1, (939, 1381)), #2630, 3071
    62 => (1, (2038, 2491)), #2038, 2491
    38 => (2, (114, 566)), #no 38b
    12 => (3, (186, 643)), #1914, 2361
    6 => (3, (790, 1246)), #1311, 1765
    21 => (4, (69, 524)), #2623, 3079
    16 => (4, (574, 1033)), #2114, 2572
    18 => (4, (1092, 1546)), #1092, 1546
    )


# we iterate on the bh sections (the dict) and we unpack the channel id (value[1]) as well as the lower limit (values[2][1]) and higher limit (values[2][1])

for (key, value) in limits_bh
    println("")
    print("This is BH:"); println(key)
    #println(value)
    println("this is channel $(value[1])")
    println("BH starts at $(value[2][1])")
    println("BH ends at $(value[2][2])")

    mT_, vt_, vL_ = get_temps_wide(
        channel = value[1],
        min_length = value[2][1], #example: min=1000 max=2000 => middle=(1000+2000)/2=1500
        #max_length = (value[2][1] + value[2][2])/2, #to the middle point
        max_length = value[2][2], #to the middle point
        min_timestamp = min_date_analysis,
        max_timestamp = max_date_analysis
    )
    
    n_rows, n_columns = size(mT_)
    default_R = mean(values(variance_dict)) #LN: This makes default_R an acceptable backup value.
    R = [get(variance_dict, string(l), default_R) for l in vL_] #LN: This uses default_R in case of missing lengths in the variance dict (e.g. 703,77 and 703,78 exists, but the one its looking for might be 703,777)
    measurement_noise = extend_variance(R, n_rows) # Extend the column variances across rows

    mT_kalman = kalman_filter(mT_, process_noise, measurement_noise, A) #LN: Kalman filtered Temperature Matrix
    
    gr() # static figures, good for report
    #pl_heatmap = heatmap(vt_, vL_, mT_kalman', title = "Borehole number $(key)", size=(1240, 900), margin=10mm, xlabel = "Time (days)", ylabel = "DTS-length within Borehole (m)", clabel = "Temperature (°C)"); savefig(pl_heatmap, "C:\\Users\\zezro\\Desktop\\MJ146X\\Plots\\HeatmapBH$(key).png")

    #averagiong wrt. time
    vT_av = mean(mT_kalman, dims=2) #LN: vT averaged with respect to length.
    vT_av_t = vec(mean(mT_kalman, dims=1)) #LN: vT averaged with respect to time.
    pl_av_T = plot(vt_, vT_av, title = "Borehole number $(key)", size=(1240, 900), margin=10mm, xlabel = "Time (days)", ylabel = "Temperature (°C)", legend = false); savefig(pl_av_T, "C:\\Users\\zezro\\Desktop\\MJ146X\\Plots\\avgT_BH$(key).png")
    #display(pl_av_T)
  
    bh_depth = mirror_bh(vL_, vL_[1], vL_[end])
    #println(length(bh_depth))
    #println(length(vL_))
    #println(length(vT_av))
    #println(length(vt_))
    
    #depthTempPlot = plot(vT_av_t, bh_depth, title = "DepthTemp borehole number $(key)", size=(1240, 900), margin=10mm, xlabel = "Temperature (°C)", ylabel = "Depth (m)", legend = false); savefig(depthTempPlot, "C:\\Users\\zezro\\Desktop\\MJ146X\\Plots\\depthPlot_average_allTimeBH$(key).png")

end
#t
#TDL: Continue with teh identification, populate the limits_bh by hand
# easiest way to tackle the "go and return" is to pick the middle point of the lower and higher bound and just trim from lower bound to middle

#plot! -> superimpose the different BH on top of each other
sp105 =superpos_plot(channel=1, min_length=174, max_length=626); #BH105
savefig(sp105, "C:\\Users\\zezro\\Desktop\\MJ146X\\Plots\\Superpositioned_Plot_BH105.png")
sp38 =superpos_plot(channel=2, min_length=114, max_length=566); #BH38
savefig(sp38, "C:\\Users\\zezro\\Desktop\\MJ146X\\Plots\\Superpositioned_Plot_BH38.png")
sp16 =superpos_plot(channel=4, min_length=574, max_length=1033) #BH16
savefig(sp16, "C:\\Users\\zezro\\Desktop\\MJ146X\\Plots\\Superpositioned_Plot_BH16.png")

#bortkommenterat:
#heatmap, plot, scatterplot