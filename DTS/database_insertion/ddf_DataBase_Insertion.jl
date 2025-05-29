using DataFrames
using LibPQ
using IniFile
using Dates
using Random
using BenchmarkTools

# ---------------------------------------------------------------------
# Automatically convert a string to an appropriate type.
function auto_convert(s::String)
    s = strip(s) # cleans up the string
    if isempty(s) #if the string is empty
        return s
    end
    # Try to parse as an Int.
    try
        return parse(Int, s)
    catch #IF To avoid the program crashing
    end
    # Try to parse as a Float (replace comma with dot).
    try
        return parse(Float64, replace(s, "," => "."))
    catch
    end
    # Try to parse as a DateTime in the "yyyy/mm/dd HH:MM:SS" format.
    try
        return DateTime(s, "yyyy/mm/dd HH:MM:SS")
    catch
    end
    # Try to parse as a Date in the "yyyy/mm/dd" format.
    try
        return Date(s, "yyyy/mm/dd")
    catch
    end
    # Try to parse as a Time in the "HH:MM:SS" format.
    try
        return Time(s, "HH:MM:SS")
    catch
    end
    return s
end

# ---------------------------------------------------------------------
# Function to parse metadata from a vector of lines.
# Always splits on tab.
function parse_meta_data(lines::Vector{<:AbstractString})
    meta = Dict{String, String}()
    for line in lines
        line = strip(line)
        if isempty(line)
            continue
        end
        # Always split on tab.
        parts = split(line, "\t")
        if length(parts) < 2
            continue
        end
        # Remove a trailing colon from the key (if present)
        key = strip(replace(parts[1], r":$" => ""))
        key = replace(key, "(" => "", ")" => "", "/" => "")
        value = strip(parts[2])
        meta[key] = value
    end
    return meta
end

# ---------------------------------------------------------------------
# Main function to parse the .ddf file.
function parse_ddf_file(filepath::String)
    # Read the whole file as a string and remove any degree symbols.
    contents = read(filepath, String)
    # Remove degree symbol in various encodings.
    contents = replace(contents, "\u00B0" => "")
    contents = replace(contents, "\xb0" => "")
    
    # Split the file into lines (handles various newline formats).
    lines = split(contents, r"[\r\n]+")
    
    # Find the header line containing "length (m)".
    header_idx = findfirst(x -> occursin("length (m)", x), lines)
    if header_idx === nothing
        error("Header not found in file.")
    end
    
    # Lines before the header are metadata.
    meta_lines = lines[1:header_idx-1]
    # Lines after the header are the time series data.
    data_lines = lines[(header_idx+1):end]
    
    # Parse the metadata into a dictionary.
    meta = parse_meta_data(meta_lines)
    
    # --- Extract and parse channel and timestamp from metadata ---
    if !haskey(meta, "forward channel")
        error("Missing 'forward channel' in metadata.")
    end
    channel_str = meta["forward channel"]
    mchan = match(r"channel\s*(\d+)", channel_str)
    if mchan === nothing
        error("Cannot extract channel number from: $channel_str")
    end
    channel = parse(Int, mchan.captures[1])
    
    if !haskey(meta, "date") || !haskey(meta, "time")
        error("Missing 'date' or 'time' in metadata.")
    end
    dt_str = string(meta["date"], " ", meta["time"])
    timestamp = DateTime(dt_str, "yyyy/mm/dd HH:MM:SS")
    
    # Remove the unparsed keys and add the parsed ones.
    delete!(meta, "forward channel")
    delete!(meta, "date")
    delete!(meta, "time")
    meta["channel"] = string(channel)
    meta["timestamp"] = string(timestamp)
    
    # Lowercase all metadata keys, remove any unit substrings (e.g. " (C)"),
    # replace spaces with underscores, and auto-convert the values.
    newmeta = Dict{Symbol, Any}()
    for (k, v) in meta
        # Remove any substring like " (C)" that might remain.
        key_no_units = replace(k, r"\s*\(\s*C\s*\)" => "")
        new_key = Symbol(lowercase(replace(key_no_units, " " => "_")))
        newmeta[new_key] = auto_convert(v)
    end
    meta_df = DataFrame(newmeta)
    
    # --- Parse the time series data ---
    # Prepare a vector to hold tuples:
    # (channel, timestamp, length, temperature, stokes, anti_stokes)
    data = Vector{Tuple{Int, DateTime, Float64, Float64, Float64, Float64}}()
    for line in data_lines
        values = split(line, r"\s+")
        if length(values) == 4
            try
                # Replace comma with dot for proper floating-point parsing.
                len_val     = parse(Float64, replace(values[1], ',' => '.'))
                temp        = parse(Float64, replace(values[2], ',' => '.'))
                stokes      = parse(Float64, replace(values[3], ',' => '.'))
                anti_stokes = parse(Float64, replace(values[4], ',' => '.'))
                
                # Filter based on channel-specific length ranges.
                if (channel == 1 && -45.100 <= len_val <= 3226.958) || #LN: Correct lengths for each channel based on value stability measured by eye.
                   (channel == 2 && -46.115 <= len_val <= 636.707)  ||
                   (channel == 3 && -46.115 <= len_val <= 2550.233) ||
                   (channel == 4 && -47.130 <= len_val <= 3139.713)
                    push!(data, (channel, timestamp, len_val, temp, stokes, anti_stokes))
                end
            catch
                continue
            end
        end
    end
    
    data_df = DataFrame(data, [:channel, :timestamp, :length, :temperature, :stokes, :anti_stokes])
    return data_df, meta_df
end


#LN: OLD: parse_ddf_file("data\\channel 1 20170801 051822 00001.ddf")

#Lägger in alla rader i Julia i en PostgreSQL-tabell
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
    
    # Debug prints:
    # println("SQL Query:")
    # println(sql_insert)
    # println("Values:")
    # println(values)
    
    try
        LibPQ.execute(conn, "BEGIN;")  # Start transaction
        LibPQ.execute(conn, sql_insert, values)
        LibPQ.execute(conn, "COMMIT;")  # Commit transaction
    catch e
        LibPQ.execute(conn, "ROLLBACK;")  # Rollback on error
        println("Error inserting data into $table: ", e)
    end
end



# Get all .ddf files
function get_all_ddf_files(folder_path)
    ddf_files = String[]
    ddf_folders = Set{String}()  #LN: Use a Set to store unique folder paths
    for (root, dirs, files) in walkdir(folder_path)
        for file in files
            if endswith(file, ".ddf")
                push!(ddf_files, joinpath(root, file))
                push!(ddf_folders, root)  #LN: Add the folder path to the set
            end
        end
    end

    #LN: Print all folders containing .ddf files
    println("Folders containing .ddf files:")
    for folder in ddf_folders
        println(folder)
        flush(stdout)
    end

    #LN: Print all file paths found
    println("All .ddf file paths:")
    for file in ddf_files
        println(file)
        flush(stdout)
    end

    return ddf_files
end



# Database connection setup
filename = "config.ini"
ini = read(Inifile(), filename)
username = get(ini, "Database", "username", "default_value")
password = get(ini, "Database", "password", "default_value")
host = get(ini, "Database", "host", "default_value")
port = get(ini, "Database", "port", "default_value")
database = get(ini, "Database", "database", "default_value")
sslmode = get(ini, "Database", "sslmode", "default_value")

conn_string = "host=$host port=$port dbname=$database user=$username password=$password sslmode=$sslmode"
conn = LibPQ.Connection(conn_string)


# Define path of the source folder and table in the PostgreSQL database
# folder_path = "C:\\Users\\zezro\\Desktop\\MJ146X\\TemperatureByMeterFibre\\Frescati_0136" #### TO BE CHANGED WITH THE PATH TO TO THE REAL FOLDER
folder_path = "C:\\Users\\zezro\\Desktop\\MJ146X\\TemperatureByMeterFibre\\Frescati_0136"
                     #LN: ^Is this AV or not...
                     #LN: Changing "data" to the real folder is difficult. 
                     #    For it to be in the same directory as the data one it needs to be in the same as the github right?

println("Starting to select all ddf files stored in the folder")
all_ddf_files = get_all_ddf_files(folder_path)
selected_files = all_ddf_files

n_files = length(selected_files)

table_name1 = "frescati_data"
table_name2 = "frescati_metadata"

# execute(conn, """
# DROP TABLE IF EXISTS $table_name1;
# """)

execute(conn, """
CREATE TABLE IF NOT EXISTS $table_name1 (
    channel INT,
    timestamp TIMESTAMP,
    length FLOAT,
    temperature FLOAT,
    stokes FLOAT,
    anti_Stokes FLOAT
    );
""")

result = execute(conn, """SELECT column_name, data_type
FROM information_schema.columns
WHERE table_schema = 'public'
AND table_name   = '$table_name1';""")
cols_data = DataFrame(result)



execute(conn, """
CREATE TABLE IF NOT EXISTS frescati_metadata (
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



# Test with batch size of 10 and measure time taken
batch_size = 2 # AV: too much data was makign it fail
# LN: What makes the time aspect important here? Isn't it unneccessary to cut 1/10 of the processing?
# AV: the 10 is in order not to fill up the RAM of the computer entirely, we process just 10 files at a time I believe
println("Starting to parse each file and create dataframe")
@time begin
    data_frames1 = DataFrame[]  # Store parsed DataFrames
    data_frames2 = DataFrame[]
    for i in 1:batch_size:length(selected_files)
        println(round(Int, 100i/length(selected_files)))
        batch_files = selected_files[i:min(i + batch_size - 1, end)]
        
        # Parse each file in the batch
        for file in batch_files
            data_df, metadata_df = parse_ddf_file(file)
            if data_df !== nothing
                push!(data_frames1, data_df)
            end
            # rename!(metadata_df, Dict(
            #     "DTS Sentinel unit serial number:" => "DTS_Sentinel_unit_serial_number",
            #     "Hardware model number:" => "Hardware_model_number",
            #     "Multiplexer serial number:" => "Multiplexer_serial_number",
            #     "Software version number:" => "Software_version_number",
            #     "T ext. ref 1 (\xb0C)" => "T_ext_ref_1_C",
            #     "T ext. ref 2 (\xb0C)" => "T_ext_ref_2_C",
            #     "T internal ref (\xb0C)" => "T_internal_ref_C",
            #     "data status" => "data_status",
            #     "default loss term (dB/km)" => "default_loss_term_dB_per_km",
            #     "differential loss correction" => "differential_loss_correction",
            #     "fibre end" => "fibre_end",
            #     "forward acquisition time" => "forward_acquisition_time",
            #     "forward channel" => "forward_channel",
            #     "k external" => "k_external",
            #     "k internal" => "k_internal",
            #     "temperature offset calibration" => "temperature_offset_calibration",
            #     "temperature slope calibration" => "temperature_slope_calibration",
            #     "reverse channel" => "reverse_channel",
            #     "multiplexer offset coefficient" => "multiplexer_offset_coefficient",
            #     "multiplexer slope coefficient" => "multiplexer_slope_coefficient",
            #     "reverse acquisition time" => "reverse_acquisition_time"
            # ))

            #if metadata_df !== nothing
            #    push!(data_frames2, metadata_df)
            #end
        end
        # Insert parsed batch data
        if !isempty(data_frames1)
            batch_data1 = vcat(data_frames1...)
            # display(batch_data1[1:10, :])
            insert_dataframe_to_postgresql(conn, batch_data1, table_name1)
            empty!(data_frames1)
        end

        # # Insert parsed batch data
        # if !isempty(data_frames2)
        #     batch_data2 = vcat(data_frames2...)
        #     insert_dataframe_to_postgresql(conn, batch_data2, table_name2)
        #     print(".")
        #     empty!(data_frames2)
        # end
    end
end

# Close the connection
close(conn)

