using DataFrames
using CSV
using Plots
using Statistics

filepath=String[]
for i in 1:3
    push!(filepath, "C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\DeviationsMeasurements\\kth $i.csv")
end
for i in 5:13
    push!(filepath, "C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\DeviationsMeasurements\\kth $i.csv")
end
# Load CSV while skipping initial lines
dfs = Vector{DataFrame}(undef, 12)
for i in 1:12
    df_raw=DataFrame[]
    df_raw= CSV.File(filepath[i]; skipto=17) |> DataFrame  # Adjust skipto accordingly
    col_names = String[]
    for col in 1:ncol(df_raw)
        push!(col_names, df_raw[1, :][col])  # Convert DataFrameRow to an array of Strings
    end
    rename!(df_raw, Symbol.(col_names))  # Convert to Symbols and set as column names
    dfs[i] = df_raw[2:end-2, :]
end

Tg=Vector{Vector{Float64}}(undef, 12)
depth=Vector{Vector{Float64}}(undef,12)
elev=Vector{Vector{Float64}}(undef,12)
for i in 1:12
    Tg[i] = parse.(Float64, dfs[i]."Temp.")
    depth[i] = parse.(Float64, dfs[i]."Depth")
    elev[i]=parse.(Float64, dfs[i]."Elev.")
end

# max reached at 8.44062e7 seconds for BH1 at time_values1[71214] and 8.15358e7 seconds for BH10 at time_values2[66427] considering the best fit sine
# min reach after 6.86382e7 seconds for BH1 at time_values1[44945] and after 6.57678e7 seconds for BH10 at time_values2[40158]
spline_BH10_max = Spline1D(fiber_lengths2[30:end-2], kalman_BH10[66427,30:end-2], k=3)
spline_BH1_max = Spline1D(fiber_lengths1[7:end-2], kalman_BH1[71214,7:end-2], k=3)
spline_BH10_min = Spline1D(fiber_lengths2[30:end-2], kalman_BH10[40158,30:end-2], k=3)
spline_BH1_min = Spline1D(fiber_lengths1[7:end-2], kalman_BH1[44945,7:end-2], k=3)
new_date = start_date + Second(time_values1[28903])

l2_ = collect(range(fiber_lengths2[30], fiber_lengths2[end-2], length = length(fiber_lengths2[30:end-2])))
L2_ = l2_ .- l2_[1]
l1_ = collect(range(fiber_lengths1[8], fiber_lengths1[end-2], length = length(fiber_lengths1[7:end-2])))
L1_ = l1_ .- l1_[1]

plotly()

plot(-elev[9], Tg[9], xlabel="Depth", size=(800,600), ylabel="Temperature", label="BH10 undisturbed")
vline!([43], linestyle=:dash, color=:black, label="Water Table")
plot!(L2_, l -> spline_BH10_max(l+fiber_lengths2[30]), label="BH10 2022-08-01", lw=2, color=:green)
plot!(L2_, l -> spline_BH10_min(l+fiber_lengths2[30]), label="BH10 2022-01-31", color=:green)
savefig("BH10_undisturbed.html")

plot(-elev[1], Tg[1], label="BH1 undisturbed", xlabel="Depth", size=(800,600), ylabel="Temperature")
plot!(L1_, l -> spline_BH1_max(l+fiber_lengths1[8]), label="BH1 2022-09-03", lw=2, color=:green)
vline!([43], linestyle=:dash, color=:black, label="Water Table")
plot!(L1_, l -> spline_BH1_min(l+fiber_lengths1[8]), label="BH1 2022-03-05", color=:green)
savefig("BH1_undisturbed.html")

Tave=Vector{Float64}(undef,12)
for j in 1:12
    for i in 1:length(Tg[j])
        if elev[j][i] < -43
            Tave[j]=mean(Tg[j][i:end])
            break
        end
    end
end
Tave[9]
Tave[1]

all_values = Float64[]
for j in 1:12
    for i in 1:length(Tg[j])
        if elev[j][i] < -43
            push!(all_values, Tg[j][i])
        end
    end
end
mean_value = mean(all_values)
