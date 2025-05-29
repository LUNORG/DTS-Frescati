### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 90a811bd-d49b-4e1a-8af7-c4c694bc35eb
using DataFrames

# Starting from 150m (a position reasonably part of the fiber in the borehole) the max length limit was found considering absolute temperature limits -15°C<T<40°C

# ╔═╡ 88af3c31-62b3-49a2-88d0-6f6f60d655ef
using LibPQ

# ╔═╡ e5baf758-1028-4f98-8a44-f4f2f34c85b7
using IniFile

# ╔═╡ 9c12469f-c462-4c04-ba15-2b585fbce373
using Dates

# ╔═╡ cb30dabe-a0ea-43dd-86f2-09e4df8d02f2
using Random

# ╔═╡ 4f191c3d-f26c-43a4-a6a9-a9ba3e742a5b
using BenchmarkTools  

# Function to parse a single .ddf file and return data and metadata DataFrames

# ╔═╡ 3fca9ef5-11e8-431b-8999-a444c8394d8c
using PlutoUI

# ╔═╡ 37179a6c-b1a5-40dd-9583-01c7b5a679b0
@bind ddf_file FilePicker()

# ╔═╡ a86482a1-8155-430c-aedd-f95b28e37b3d
readline(ddf_file["data"])

# ╔═╡ 064fce25-0d0f-4e28-9454-d741329bc501
byte_data = reinterpret(UInt8, ddf_file["data"])

# ╔═╡ 213e1d3c-5e6b-4b30-9ed2-d6c848697d65
io = IOBuffer(byte_data)

# ╔═╡ 97c2b66e-2954-4dd6-9103-f1f7ad71cef8
eachline(io)

# ╔═╡ 84048362-be71-4f60-89e1-4f8cd008b5a2
temp_file = tempname()

# ╔═╡ 647cd907-488d-478c-9114-6e1da9683a2c
open(temp_file, "r") do io_
	for line in eachline(io_)
		println(line)
	end
end

# ╔═╡ 574044f4-f175-407b-8d20-4937dc38a676
function parse_ddf_file(filepath)
    try
        lines = readlines(filepath)
        lines = replace.(lines, "°" => "")  # Remove "°" character

        # Extract metadata from the first lines
        split_lines = split(lines[1], r"[\t\r]")    #LN: ska inte bara vara 1a, utan rad 1-26...
        metadata_lines = []                         #LN: "split_lines = split(join(lines[1:26], " "), r"[\t\r]")"
        for i in 1:length(split_lines)              
            if split_lines[i] == "length (m)"
                break
            end
            push!(metadata_lines, split_lines[i])
        end
        
        # Create metadata DataFrame
        metadata_columns = metadata_lines[1:2:end]
        values = metadata_lines[2:2:end]
        metadata_df = DataFrame(Dict(metadata_columns .=> values))

        data = []
        # Determine the channel from metadata
        channel_forward_value = metadata_df[!, "forward channel"][1]
        channel = channel_forward_value == "channel 1" ? 1 :
                  channel_forward_value == "channel 2" ? 2 :
                  channel_forward_value == "channel 3" ? 3 :        #LN: ändrat från bara 2 channels
                  channel_forward_value == "channel 4" ? 4 : NaN

        # Process data lines
        for line in lines[2:end]
            values = split(line)
            if length(values) == 4  # Ensure 4 columns
                try
                    length_val = parse(Float64, values[1])
                    temperature = parse(Float64, values[2])
                    stoke = parse(Float64, values[3])
                    anti_stokes = parse(Float64, values[4])

                    #LN: TIDIGARE KOD: if (channel == 1 && -50.678 <= length_val <= 1059.289) || (channel == 2 && -1.977 <= length_val <= 511.408)
                    #LN: Ny kod för 4 channels:
                    if (channel == 1 && -45.100 <= length_val <= 3226.958) ||
                       (channel == 2 && -46.115 <= length_val <= 636.707) ||
                       (channel == 3 && -46.115 <= length_val <= 2550.233) ||
                       (channel == 4 && -47.130 <= length_val <= 3139.713)
                    #LN: slut ny kod. Kräver dock att man vet de olika sträckorna för de olika channels, 
                    #    vilket jag tog reda på genom att läsa ddt filens logg.
                        push!(data, (channel, length_val, temperature, stoke, anti_stokes))
                    end
                catch
                    continue  # Skip lines that cannot be parsed
                end
            end
        end

        # Create data DataFrame
        data_df = DataFrame(data, [:channel, :length, :temperature, :Stokes, :anti_Stokes])
        
        # Rebuild the DataFrames with datetime as the first column
        datetime_col = DateTime(metadata_df.date[1] * " " * metadata_df.time[1], dateformat"yyyy/mm/dd HH:MM:SS")
        metadata_df[!, :datetime] = datetime_col
        metadata_df = select(metadata_df, :datetime, Not([:datetime, :date, :time]))
        
        expanded_datetime_col = fill(datetime_col, nrow(data_df))
        data_df[!, :datetime] = expanded_datetime_col
        data_df = select(data_df, :datetime, Not(:datetime))

        return data_df, metadata_df
    catch e
        println("Error parsing DTS file: $e")
        return nothing, nothing
    end
    #println("Processing file: $filepath, Assigned Channel: $channel")  #LN: detta är för att testa ifall det sparas korrekt.
    #SELECT DISTINCT channel FROM DTSold_data;
end

# Function to insert all rows in a DataFrame within a single transaction

# ╔═╡ 8eb31a99-cd72-4cc2-9bdf-948ee51aeeed
function insert_dataframe_to_postgresql(conn::LibPQ.Connection, df::DataFrame, table::String)
    columns = join(names(df), ", ")
    sql_insert = "INSERT INTO $table ($columns) VALUES "
    
    num_columns = ncol(df)
    num_rows = nrow(df)
    placeholders_list = []
    for i in 1:num_rows
        row_placeholders = join(["\$" * string((i - 1) * num_columns + j) for j in 1:num_columns], ", ")
        push!(placeholders_list, "($row_placeholders)")
    end
    placeholders = join(placeholders_list, ", ")
    sql_insert *= placeholders * ";"
    
    values = collect(Iterators.flatten(eachrow(df)))

    try
        LibPQ.execute(conn, "BEGIN;")  # Start transaction
        LibPQ.execute(conn, sql_insert, values)
        LibPQ.execute(conn, "COMMIT;")  # Commit transaction
    catch e
        LibPQ.execute(conn, "ROLLBACK;")  # Rollback on error
        println("Error inserting data into $table: ", e)
    end
end

# Database connection setup

# ╔═╡ c6a6d1e7-f346-4ee0-85bc-9e232dc941f0
filename = "config.ini"

# ╔═╡ 00644533-8777-43a3-b2e7-caa361ddd11d
ini = read(Inifile(), filename)

# ╔═╡ b4a2c21e-5f17-4f7b-a528-910e66d2e2bf
username = get(ini, "Database", "username", "default_value")

# ╔═╡ 4ff220e0-4a82-4f05-95c1-f64bcb45d256
password = get(ini, "Database", "password", "default_value")

# ╔═╡ c41cfa2a-8258-47cd-8ee8-141ca5c3da46
host = get(ini, "Database", "host", "default_value")

# ╔═╡ b4de1433-c5e2-4d9b-a6ee-c36deae2d222
port = get(ini, "Database", "port", "default_value")

# ╔═╡ d0a37bc1-b896-4a5d-a054-cf7773eed1f9
database = get(ini, "Database", "database", "default_value")

# ╔═╡ 2e198c39-4756-4d24-a67f-9eb3b28c4eb6
sslmode = get(ini, "Database", "sslmode", "default_value")

# ╔═╡ 44385096-9ad1-4b12-b1ef-06e36461f8d4
conn_string = "host=$host port=$port dbname=$database user=$username password=$password sslmode=$sslmode"

# ╔═╡ 4db9ac53-a38c-4aa9-a562-6b77c74ad63b
conn = LibPQ.Connection(conn_string)

# Define path of the source folder and table in the PostgreSQL database

# ╔═╡ 0def2d7c-090d-4d8e-8bbf-1407dbca2391
folder_path = "DTS_data"

# ╔═╡ 32b1275e-ed67-486b-bf8b-5716b1fc0494
table_name1 = "DTSold_data"

# ╔═╡ 61f6ed6e-088d-4bbd-9299-5ce010a44336
table_name2 = "DTSold_metadata"

# ╔═╡ 71925f6c-8ca0-4645-9502-687cf97e4694
execute(conn, """
CREATE TABLE IF NOT EXISTS DTSold_data (
    datetime TIMESTAMP,
    channel FLOAT,
    length FLOAT,
    temperature FLOAT,
    Stokes FLOAT,
    anti_Stokes FLOAT
    );
""")

# ╔═╡ 72bb690b-568c-477e-9446-3ca0213d91ae
execute(conn, """
CREATE TABLE IF NOT EXISTS DTSold_metadata (
    datetime TIMESTAMP,
    DTS_Sentinel_unit_serial_number TEXT,
    Multiplexer_serial_number TEXT,
    Hardware_model_number TEXT,
    Software_version_number TEXT,
    data_status TEXT,
    installation TEXT,
    differential_loss_correction TEXT,
    forward_channel TEXT,
    reverse_channel TEXT,
    forward_acquisition_time FLOAT,
    reverse_acquisition_time FLOAT,
    gamma FLOAT,
    k_internal FLOAT,
    k_external FLOAT,
    temperature_offset_calibration FLOAT,
    default_loss_term_dB_per_km FLOAT,
    temperature_slope_calibration FLOAT,
    multiplexer_offset_coefficient FLOAT,
    multiplexer_slope_coefficient FLOAT,
    fibre_end TEXT,
    T_internal_ref_C FLOAT,
    T_ext_ref_1_C FLOAT,
    T_ext_ref_2_C FLOAT
);
""")

# Get all .ddf files

# ╔═╡ a2f906f2-fc55-46f9-aff5-6b9dc3727a45
function get_all_ddf_files(folder_path)
    ddf_files = String[]
    for (root, dirs, files) in walkdir(folder_path)
        for file in files
            if endswith(file, ".ddf")
                push!(ddf_files, joinpath(root, file))
            end
        end
    end
    return ddf_files
end

# ╔═╡ 7cc450c0-7f23-43a6-981f-39669706ef73
println("Starting to select all ddf files stored in the folder")

# ╔═╡ 5292d7ab-1ea4-4ac0-9b5d-2b29a5d17a22
all_ddf_files = get_all_ddf_files(folder_path)

# ╔═╡ 4e8272b6-c9a1-4b5e-8442-0ba48f373c53
selected_files = all_ddf_files

# Test with batch size of 10 and measure time taken

# ╔═╡ 67479a49-9948-4779-b322-a93a945544ec
batch_size = 10
# LN: What makes the time aspect important here? Isn't it unneccessary to cut 1/10 of the processing?

# ╔═╡ f5812947-2bc9-4abf-aef6-ad9a0b30e55a
println("Starting to parse each file and create dataframe")

# ╔═╡ 793c50f2-8dac-4986-90da-0bb9c6e1141e
@time begin
    data_frames1 = DataFrame[]  # Store parsed DataFrames
    data_frames2 = DataFrame[]
    for i in 1:batch_size:length(selected_files)
        batch_files = selected_files[i:min(i + batch_size - 1, end)]
        
        # Parse each file in the batch
        for file in batch_files
            data_df, metadata_df = parse_ddf_file(file)

            if data_df !== nothing
                push!(data_frames1, data_df)
            end

            rename!(metadata_df, Dict(
                "DTS Sentinel unit serial number:" => "DTS_Sentinel_unit_serial_number",
                "Hardware model number:" => "Hardware_model_number",
                "Multiplexer serial number:" => "Multiplexer_serial_number",
                "Software version number:" => "Software_version_number",
                "T ext. ref 1 (\xb0C)" => "T_ext_ref_1_C",
                "T ext. ref 2 (\xb0C)" => "T_ext_ref_2_C",
                "T internal ref (\xb0C)" => "T_internal_ref_C",
                "data status" => "data_status",
                "default loss term (dB/km)" => "default_loss_term_dB_per_km",
                "differential loss correction" => "differential_loss_correction",
                "fibre end" => "fibre_end",
                "forward acquisition time" => "forward_acquisition_time",
                "forward channel" => "forward_channel",
                "k external" => "k_external",
                "k internal" => "k_internal",
                "temperature offset calibration" => "temperature_offset_calibration",
                "temperature slope calibration" => "temperature_slope_calibration",
                "reverse channel" => "reverse_channel",
                "multiplexer offset coefficient" => "multiplexer_offset_coefficient",
                "multiplexer slope coefficient" => "multiplexer_slope_coefficient",
                "reverse acquisition time" => "reverse_acquisition_time"
            ))

            if metadata_df !== nothing
                push!(data_frames2, metadata_df)
            end
        end

        # Insert parsed batch data
        if !isempty(data_frames1)
            batch_data1 = vcat(data_frames1...)
            insert_dataframe_to_postgresql(conn, batch_data1, table_name1)
            print(".")
            empty!(data_frames1)
        end

        # Insert parsed batch data
        if !isempty(data_frames2)
            batch_data2 = vcat(data_frames2...)
            insert_dataframe_to_postgresql(conn, batch_data2, table_name2)
            print(".")
            empty!(data_frames2)
        end
    end
end

# Close the connection

# ╔═╡ 4cc5335d-6364-415e-886a-4a4ff9ee5331
close(conn)

# ╔═╡ 94d535e6-1097-4ce3-ac7f-18b668e08f49
max_length = parse_all_ddf_in_folder("DTS_data\\full data set\\KTH_LiL_BH10_12")

# ╔═╡ 2ca81abe-bfa1-4b7f-985b-b87acc3207b6
ddf_file

# ╔═╡ c9ac82a4-d5d5-41f0-87a8-f1d3732b3bb9
string_ = ddf_file["data"]

# ╔═╡ fbf2d4ad-2114-49e2-9de2-a9a71d41e8d6
content = String(ddf_file["data"])

# ╔═╡ 0836d8ca-46b6-4133-b8b9-c9dd4cf7414f
data = []
        # Determine the channel from metadata

# ╔═╡ 3ae07cf5-57f7-4f20-bafd-c3b3a86183b1
begin
	lines = split(content, r"[\r]")
	for line in lines[26:end]  # Skip the first line with metadata
            values = split(line, r"[\t\n]")[2:end]
            if length(values) == 4  # Ensure there are exactly 4 columns
                try
                    length_val = parse(Float64, values[1])
                    temperature = parse(Float64, values[2])
                    push!(data, (length_val, temperature))  # Store data
                catch parse_error
                    println("Error parsing values in data line: $line")
                    continue  # Skip lines that cannot be parsed
                end
            end
        end
        
        # Create data DataFrame
end

# ╔═╡ 37ef906b-797e-45c4-baff-c2131e3e80cd
lines[28]

# ╔═╡ 6bd5c132-7b54-4677-bfe2-b20d0693055b
v = split(lines[168], r"[\t\n]")[2:end]

# ╔═╡ 6276c9e6-4b3d-4c55-9af5-dde15f1179c6
split_lines = split(lines[26], r"[\t\r]")  # Split first line to extract metadata

# ╔═╡ c6784adb-c9b0-455b-b7c0-32444a672a15
channel_forward_value = metadata_df[!, "forward channel"][1]

# ╔═╡ 417d2e85-477c-4d64-991c-0e091807a5a4
channel = channel_forward_value == "channel 1" ? 1 :
                channel_forward_value == "channel 2" ? 2 : NaN

# ╔═╡ 5e27d535-4480-48d5-a5bf-e4b38757430c
for line in lines[2:end]  # Skip the first line with metadata
            values = split(line)
            if length(values) == 4  # Ensure there are exactly 4 columns
                try
                    length_val = parse(Float64, values[1])
                    temperature = parse(Float64, values[2])
                    push!(data, (length_val, temperature))  # Store data
                catch parse_error
                    println("Error parsing values in data line: $line")
                    continue  # Skip lines that cannot be parsed
                end
            end
        end
        
        # Create data DataFrame

# ╔═╡ 6963cd34-e912-4d86-91b2-16b262023221
data_df = DataFrame(data, [:length, :temperature])

# ╔═╡ 12b35aba-45f6-4cce-911d-27f6ad47cb0a
filtered_data_df = filter_dataframe(data_df)
        # Find length that meets the condition and update max_length if needed

# ╔═╡ 836d99f3-b241-47b1-badf-13cdc43f9d1a
length_value = find_length_for_temperature_jump(filtered_data_df)

# ╔═╡ 8700b513-6d1c-405d-b442-56f53fbf073e
if channel == 1 && length_value != nothing && (isnothing(max_length_channel_1[]) || length_value > max_length_channel_1[])
            max_length_channel_1[] = length_value  # Update the max_length if a new one is found
        elseif channel == 2 && length_value != nothing && (isnothing(max_length_channel_2[]) || length_value > max_length_channel_2[])
            max_length_channel_2[] = length_value
        end

# ╔═╡ 9016055f-e236-4258-be8d-083b60d774c1
function parse_ddf_file_pluto(file_pointer, max_length_channel_1=nothing, max_length_channel_2=nothing)
    try
		#AV: NEW STUFF
        #lines = readlines(filepath)
		content = String(file_pointer["data"])
		lines = split(content, r"[\r]")
		# ---
        split_lines = split(lines[1], r"[\t\r]")  # Split first line to extract metadata
        metadata_lines = []
        
        for i in 1:length(split_lines)
            if split_lines[i] == "length (m)"
                break  # Stop when "length (m)" is found
            end
            push!(metadata_lines, split_lines[i])
        end
        
        # Create metadata DataFrame
        metadata_columns = metadata_lines[1:2:end]
        values = metadata_lines[2:2:end]
        metadata_df = DataFrame(Dict(metadata_columns .=> values))
        
        # Process data lines
        data = []
        # Determine the channel from metadata
        channel_forward_value = metadata_df[!, "forward channel"][1]
        channel = channel_forward_value == "channel 1" ? 1 :
                channel_forward_value == "channel 2" ? 2 : NaN
        
        for line in lines[2:end]  # Skip the first line with metadata
            values = split(line)
            if length(values) == 4  # Ensure there are exactly 4 columns
                try
                    length_val = parse(Float64, values[1])
                    temperature = parse(Float64, values[2])
                    push!(data, (length_val, temperature))  # Store data
                catch parse_error
                    println("Error parsing values in data line: $line")
                    continue  # Skip lines that cannot be parsed
                end
            end
        end
        
        # Create data DataFrame
        data_df = DataFrame(data, [:length, :temperature])
        filtered_data_df = filter_dataframe(data_df)
        # Find length that meets the condition and update max_length if needed
        length_value = find_length_for_temperature_jump(filtered_data_df)
        if channel == 1 && length_value != nothing && (isnothing(max_length_channel_1[]) || length_value > max_length_channel_1[])
            max_length_channel_1[] = length_value  # Update the max_length if a new one is found
        elseif channel == 2 && length_value != nothing && (isnothing(max_length_channel_2[]) || length_value > max_length_channel_2[])
            max_length_channel_2[] = length_value
        end
    catch e
        println("Error parsing DTS file: $e")
        return nothing  # Return nothing if any error occurs
    end
end

# Function to crawl through a folder and parse all .ddf files

# ╔═╡ 3d2c62df-33db-4a58-86e8-d95b4e806c2e


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
IniFile = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
LibPQ = "194296ae-ab2e-5f79-8cd4-7183a0a5a0d1"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
BenchmarkTools = "~1.6.0"
DataFrames = "~1.7.0"
IniFile = "~0.5.1"
LibPQ = "~1.18.0"
PlutoUI = "~0.7.61"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.3"
manifest_format = "2.0"
project_hash = "86ef91cd99138fc3078b00d184c32a4895f9969f"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BenchmarkTools]]
deps = ["Compat", "JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "e38fbc49a620f5d0b660d7f543db1009fe0f8336"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.6.0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DBInterface]]
git-tree-sha1 = "a444404b3f94deaa43ca2a58e18153a82695282b"
uuid = "a10d1c49-ce27-4219-8d33-6db1a4562965"
version = "2.6.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fb61b4812c49343d7ef0b533ba982c46021938a6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.7.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Decimals]]
git-tree-sha1 = "e98abef36d02a0ec385d68cd7dadbce9b28cbd88"
uuid = "abce61dc-4473-55a0-ba07-351d65e31d42"
version = "0.4.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.ICU_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "20b6765a3016e1fca0c9c93c80d50061b94218b7"
uuid = "a51ab1cf-af8e-5615-a023-bc2c838bba6b"
version = "69.1.0+0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.Infinity]]
deps = ["Dates", "Random", "Requires"]
git-tree-sha1 = "cf8234411cbeb98676c173f930951ea29dca3b23"
uuid = "a303e19e-6eb4-11e9-3b09-cd9505f79100"
version = "0.2.4"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InlineStrings]]
git-tree-sha1 = "45521d31238e87ee9f9732561bfee12d4eebd52d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.2"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "ac0aaa807ed5eaf13f67afe188ebc07e828ff640"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.10.0"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.Kerberos_krb5_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0f2899fdadaab4b8f57db558ba21bdb4fb52f1f0"
uuid = "b39eb1a6-c29a-53d7-8c32-632cd16f18da"
version = "1.21.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LayerDicts]]
git-tree-sha1 = "6087ad3521d6278ebe5c27ae55e7bbb15ca312cb"
uuid = "6f188dcb-512c-564b-bc01-e0f76e72f166"
version = "1.0.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibPQ]]
deps = ["CEnum", "DBInterface", "Dates", "Decimals", "DocStringExtensions", "FileWatching", "Infinity", "Intervals", "IterTools", "LayerDicts", "LibPQ_jll", "Libdl", "Memento", "OffsetArrays", "SQLStrings", "Tables", "TimeZones", "UTCDateTimes"]
git-tree-sha1 = "3d227cd13cbf1e9a54d7748dab33e078da6f9168"
uuid = "194296ae-ab2e-5f79-8cd4-7183a0a5a0d1"
version = "1.18.0"

[[deps.LibPQ_jll]]
deps = ["Artifacts", "ICU_jll", "JLLWrappers", "Kerberos_krb5_jll", "Libdl", "OpenSSL_jll", "Zstd_jll"]
git-tree-sha1 = "09163f837936c8cc44f4691cb41d805eb1769642"
uuid = "08be9ffa-1c94-5ee5-a977-46a84ec9b350"
version = "16.0.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "1833212fd6f580c20d4291da9c1b4e8a655b128e"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.0.0"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Memento]]
deps = ["Dates", "Distributed", "Requires", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "bb2e8f4d9f400f6e90d57b34860f6abdc51398e5"
uuid = "f28f55f0-a522-5efc-85c2-fe41dfb9b2d9"
version = "1.4.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "2c140d60d7cb82badf06d8783800d0bcd1a7daa2"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.8.1"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "5e1897147d1ff8d98883cda2be2187dcf57d8f0c"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.15.0"

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

    [deps.OffsetArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a9697f1d06cc3eb3fb3ad49cc67f2cfabaac31ea"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.16+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "7e71a55b87222942f0f9337be62e26b1f103d3e4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.61"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Profile]]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SQLStrings]]
git-tree-sha1 = "55de0530689832b1d3d43491ee6b67bd54d3323c"
uuid = "af517c2e-c243-48fa-aab8-efac3db270f5"
version = "0.1.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "725421ae8e530ec29bcbdddbe91ff8053421d023"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TZJData]]
deps = ["Artifacts"]
git-tree-sha1 = "7def47e953a91cdcebd08fbe76d69d2715499a7d"
uuid = "dc5dba14-91b3-4cab-a142-028a31da12f7"
version = "1.4.0+2025a"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TimeZones]]
deps = ["Artifacts", "Dates", "Downloads", "InlineStrings", "Mocking", "Printf", "Scratch", "TZJData", "Unicode", "p7zip_jll"]
git-tree-sha1 = "38bb1023fb94bfbaf2a29e1e0de4bbba6fe0bf6d"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.21.2"
weakdeps = ["RecipesBase"]

    [deps.TimeZones.extensions]
    TimeZonesRecipesBaseExt = "RecipesBase"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UTCDateTimes]]
deps = ["Dates", "TimeZones"]
git-tree-sha1 = "4af3552bf0cf4a071bf3d14bd20023ea70f31b62"
uuid = "0f7cfa37-7abf-4834-b969-a8aa512401c2"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "622cf78670d067c738667aaa96c553430b65e269"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╠═90a811bd-d49b-4e1a-8af7-c4c694bc35eb
# ╠═88af3c31-62b3-49a2-88d0-6f6f60d655ef
# ╠═e5baf758-1028-4f98-8a44-f4f2f34c85b7
# ╠═9c12469f-c462-4c04-ba15-2b585fbce373
# ╠═cb30dabe-a0ea-43dd-86f2-09e4df8d02f2
# ╠═4f191c3d-f26c-43a4-a6a9-a9ba3e742a5b
# ╠═a86482a1-8155-430c-aedd-f95b28e37b3d
# ╠═37179a6c-b1a5-40dd-9583-01c7b5a679b0
# ╠═064fce25-0d0f-4e28-9454-d741329bc501
# ╠═213e1d3c-5e6b-4b30-9ed2-d6c848697d65
# ╠═97c2b66e-2954-4dd6-9103-f1f7ad71cef8
# ╠═84048362-be71-4f60-89e1-4f8cd008b5a2
# ╠═647cd907-488d-478c-9114-6e1da9683a2c
# ╠═574044f4-f175-407b-8d20-4937dc38a676
# ╠═8eb31a99-cd72-4cc2-9bdf-948ee51aeeed
# ╠═c6a6d1e7-f346-4ee0-85bc-9e232dc941f0
# ╠═00644533-8777-43a3-b2e7-caa361ddd11d
# ╠═b4a2c21e-5f17-4f7b-a528-910e66d2e2bf
# ╠═4ff220e0-4a82-4f05-95c1-f64bcb45d256
# ╠═c41cfa2a-8258-47cd-8ee8-141ca5c3da46
# ╠═b4de1433-c5e2-4d9b-a6ee-c36deae2d222
# ╠═d0a37bc1-b896-4a5d-a054-cf7773eed1f9
# ╠═2e198c39-4756-4d24-a67f-9eb3b28c4eb6
# ╠═44385096-9ad1-4b12-b1ef-06e36461f8d4
# ╠═4db9ac53-a38c-4aa9-a562-6b77c74ad63b
# ╠═0def2d7c-090d-4d8e-8bbf-1407dbca2391
# ╠═32b1275e-ed67-486b-bf8b-5716b1fc0494
# ╠═61f6ed6e-088d-4bbd-9299-5ce010a44336
# ╠═71925f6c-8ca0-4645-9502-687cf97e4694
# ╠═72bb690b-568c-477e-9446-3ca0213d91ae
# ╠═a2f906f2-fc55-46f9-aff5-6b9dc3727a45
# ╠═7cc450c0-7f23-43a6-981f-39669706ef73
# ╠═5292d7ab-1ea4-4ac0-9b5d-2b29a5d17a22
# ╠═4e8272b6-c9a1-4b5e-8442-0ba48f373c53
# ╠═67479a49-9948-4779-b322-a93a945544ec
# ╠═f5812947-2bc9-4abf-aef6-ad9a0b30e55a
# ╠═793c50f2-8dac-4986-90da-0bb9c6e1141e
# ╠═4cc5335d-6364-415e-886a-4a4ff9ee5331
# ╠═94d535e6-1097-4ce3-ac7f-18b668e08f49
# ╠═3fca9ef5-11e8-431b-8999-a444c8394d8c
# ╠═2ca81abe-bfa1-4b7f-985b-b87acc3207b6
# ╠═c9ac82a4-d5d5-41f0-87a8-f1d3732b3bb9
# ╠═fbf2d4ad-2114-49e2-9de2-a9a71d41e8d6
# ╠═37ef906b-797e-45c4-baff-c2131e3e80cd
# ╠═6bd5c132-7b54-4677-bfe2-b20d0693055b
# ╠═3ae07cf5-57f7-4f20-bafd-c3b3a86183b1
# ╠═6276c9e6-4b3d-4c55-9af5-dde15f1179c6
# ╠═0836d8ca-46b6-4133-b8b9-c9dd4cf7414f
# ╠═c6784adb-c9b0-455b-b7c0-32444a672a15
# ╠═417d2e85-477c-4d64-991c-0e091807a5a4
# ╠═5e27d535-4480-48d5-a5bf-e4b38757430c
# ╠═6963cd34-e912-4d86-91b2-16b262023221
# ╠═12b35aba-45f6-4cce-911d-27f6ad47cb0a
# ╠═836d99f3-b241-47b1-badf-13cdc43f9d1a
# ╠═8700b513-6d1c-405d-b442-56f53fbf073e
# ╠═9016055f-e236-4258-be8d-083b60d774c1
# ╠═3d2c62df-33db-4a58-86e8-d95b4e806c2e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
