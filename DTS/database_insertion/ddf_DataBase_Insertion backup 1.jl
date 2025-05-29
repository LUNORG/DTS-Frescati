using DataFrames
using LibPQ
using IniFile
using Dates
using Random
using BenchmarkTools  

# Function to parse a single .ddf file and return data and metadata DataFrames
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
folder_path = "DTS_data"
table_name1 = "DTSold_data"
table_name2 = "DTSold_metadata"

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

println("Starting to select all ddf files stored in the folder")
all_ddf_files = get_all_ddf_files(folder_path)
selected_files = all_ddf_files

# Test with batch size of 10 and measure time taken
batch_size = 10
# LN: What makes the time aspect important here? Isn't it unneccessary to cut 1/10 of the processing?

println("Starting to parse each file and create dataframe")
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
close(conn)