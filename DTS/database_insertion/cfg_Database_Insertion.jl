using DataFrames
using LibPQ
using IniFile
using Dates
using Random
using BenchmarkTools  

# Function to parse a single .cfg file and return data and metadata DataFrames
function parse_cfg_file(filepath)
    try
        lines = readlines(filepath)                 #LN: Vad är en cfg fil i detta sammanhang; vi har ju ingen sådan...
        split_lines = split(lines[1], r"[\t\r]")    #LN: ska inte bara vara 1a, utan rad 1-26...
                                                    #LN: "split_lines = split(join(lines[1:26], " "), r"[\t\r]")"
        metadata_lines = [] 
        found_index = nothing       
        for i in 1:length(split_lines)
            if split_lines[i] == "CHANNEL: "
                found_index = i 
                break  # Stop when "length (m)" is found
            end
            push!(metadata_lines, split_lines[i])
        end

        data_lines1 = []
        data_lines2 = []
        new_found_index = nothing

        for i in found_index+1:length(split_lines)
            if split_lines[i] == "CHANNEL: "        #LN: Söker denna efter "CHANNEL:" i datan?
                new_found_index = i                 #LN: Det ser ut så, isf ska detta ändras till "channel " för oss.
                break                               #LN: detta kan vara fel; eftersom vi har 4 mappar med nmanen channel 1-4 och de har sedan all data som egna filer inom sig.
            end                                     #LN: Läser denna ens flera filer, ser ut som att hon bara har en fil med alla data som den läser medan vi har fyra mappar med submappar och sjukt många filer inom sig var.
        end

        for i in found_index:length(split_lines)
            if i<new_found_index
                push!(data_lines1, split_lines[i])
            else
                push!(data_lines2, split_lines[i])
            end
        end
        pop!(data_lines2)

        metadata_columns = metadata_lines[1:2:end-1]
        values = metadata_lines[2:2:end]
        metadata_df = DataFrame(Dict(metadata_columns .=> values))

        data1_lines = data_lines1[1:2:end]
        values_data1 = data_lines1[2:2:end]
        data1_df = DataFrame(data1_lines .=> values_data1, makeunique=true)
        data1_combined = hcat(metadata_df, data1_df)
        select!(data1_combined, :"CHANNEL: ", Not([:"CHANNEL: "]))

        data2_lines = data_lines2[1:2:end]
        values_data2 = data_lines2[2:2:end]
        data2_df = DataFrame(data2_lines .=> values_data2, makeunique=true)
        data2_combined = hcat(metadata_df, data2_df)
        select!(data2_combined, :"CHANNEL: ", Not([:"CHANNEL: "]))

        final_df = vcat(data1_combined, data2_combined)

        # Extract filename from filepath, isolate the date string
        filename = basename(filepath)        
        insertcols!(final_df, 1, :filename => fill(filename, nrow(final_df)))
        
        return final_df
    catch e
        println("Error parsing DTS file: $e")
        return nothing
    end
end

# Folder "full data set\Tw first measurement.cfg" is the only one with different content structure so I skipped it and inserted it separately

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

function get_all_cfg_files(folder_path)
    cfg_files = String[]
    for (root, dirs, files) in walkdir(folder_path)
        for file in files
            if endswith(file, ".cfg")
                file_path = joinpath(root, file)
                push!(cfg_files, file_path)
            end
        end
    end
    return cfg_files
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

# Define path and table
folder_path = "database_inserting\\DTS_data\\full data set"
table_name = "DTSold_configuration"

execute(conn, """
CREATE TABLE IF NOT EXISTS DTSold_configuration (
    filename TEXT,
    channel TEXT,
    "dts_sentinel_unit_serial_number" TEXT,
    "hardware_model_number" TEXT,
    "multiplexer_serial_number" TEXT,
    "software_version_number" TEXT,
    "assumed_t_internal" TEXT,
    comments TEXT,
    continuous TEXT,
    "default_nga_plus_ngr_div2" TEXT,
    "default_ngs_plus_ngr_div2" TEXT,
    "default_diff_loss_term" TEXT,
    gamma TEXT,
    installation TEXT,
    "internal_lead_length" TEXT,
    "internal_reference_start" TEXT,
    "internal_reference_stop" TEXT,    
    "number_of_channels" TEXT,
    "number_of_repetitions" TEXT,
    "offset_reference_start" TEXT,
    "offset_reference_stop" TEXT,
    "repetition_time_s" TEXT,
    "save_data" TEXT,
    "channel_active" TEXT, 
    "range_in_points" TEXT,
    "sample_length" TEXT,
    "moving_average_in_points" TEXT,
    "stokes_length_correction" TEXT,
    range TEXT,
    "temperature_offset_calculation" TEXT,
    "fixed_t_offset_calibration_value" TEXT,
    "offset_reference_type" TEXT,
    "offset_reference_start_m" TEXT,
    "offset_reference_stop_m" TEXT,
    "offset_temperature_source" TEXT,
    "fixed_temperature_for_offset" TEXT,
    "temperature_slope_calculation" TEXT,
    "fixed_t_slope_calibration_value" TEXT,
    "temperature_slope_calculation_1" TEXT,
    "first_slope_reference_start" TEXT,
    "first_slope_reference_stop" TEXT,
    "first_slope_ref_temperature_source" TEXT,
    "first_slope_ref_fixed_temperature" TEXT,
    "second_slope_reference_type" TEXT,
    "second_slope_reference_start" TEXT,
    "second_slope_reference_stop" TEXT,
    "second_slope_ref_temperature_source" TEXT,
    "second_slope_ref_fixed_temperature" TEXT,
    "differential_loss_correction" TEXT,
    "paired_channel" TEXT,
    direction TEXT,
    "fibre_end_point" TEXT,
    "spatial_averaging_m" TEXT,
    "measurement_time_s" TEXT,
    "number_of_zones" TEXT
    );
""")

function rename_columns(metadataframe::DataFrame)
    rename!(metadataframe,
    "DTS Sentinel unit serial number:" => "DTS_Sentinel_unit_serial_number",
    "Multiplexer serial number:" => "Multiplexer_serial_number",
    "Hardware model number:" => "Hardware_model_number",
    "Software version number:" => "Software_version_number",
    "CHANNEL: " => "channel",   
    "Stokes length correction" => "Stokes_length_correction",
    "internal lead length" => "internal_lead_length",
    "default diff loss term" => "default_diff_loss_term",
    "default (ngA + ngR)/2" => "default_ngA_plus_ngR_div2",
    "default (ngS + ngR)/2" => "default_ngS_plus_ngR_div2",
    "offset reference start" => "offset_reference_start",
    "offset reference stop" => "offset_reference_stop",
    "offset reference start (m)" => "offset_reference_start_m",
    "offset reference stop (m)" => "offset_reference_stop_m",
    "internal reference start" => "internal_reference_start",
    "internal reference stop" => "internal_reference_stop",
    "assumed T internal" => "assumed_T_internal",
    "channel active" => "channel_active", 
    "differential loss correction" => "differential_loss_correction",
    "fibre end point" => "fibre_end_point",
    "first slope ref fixed temperature" => "first_slope_ref_fixed_temperature",
    "first slope ref temperature source" => "first_slope_ref_temperature_source",
    "first slope reference start" => "first_slope_reference_start",
    "first slope reference stop" => "first_slope_reference_stop",
    "fixed T offset calibration value" => "fixed_T_offset_calibration_value",
    "fixed T slope calibration value" => "fixed_T_slope_calibration_value",
    "fixed temperature for offset" => "fixed_temperature_for_offset",
    "save data" => "save_data",
    "repetition time (s)" => "repetition_time_s",
    "measurement time (s)" => "measurement_time_s",
    "moving average in points" => "moving_average_in_points",
    "number of channels" => "number_of_channels",
    "number of repetitions" => "number_of_repetitions",
    "number of zones" => "number_of_zones",
    "offset reference type" => "offset_reference_type",
    "offset temperature source" => "offset_temperature_source",
    "paired channel" => "paired_channel",
    "range in points" => "range_in_points",
    "sample length" => "sample_length",
    "second slope ref fixed temperature" => "second_slope_ref_fixed_temperature",
    "second slope ref temperature source" => "second_slope_ref_temperature_source",
    "second slope reference start" => "second_slope_reference_start",
    "second slope reference stop" => "second_slope_reference_stop",
    "second slope reference type" => "second_slope_reference_type",
    "spatial averaging (m)" => "spatial_averaging_m",
    "temperature offset calculation" => "temperature_offset_calculation",
    "temperature slope calculation" => "temperature_slope_calculation",
    "temperature slope calculation_1" => "temperature_slope_calculation_1" 
    )
end

println("Starting to select al ddf files stored in the folder")
all_cfg_files = get_all_cfg_files(folder_path)
selected_files = all_cfg_files

# Test with batch size of 10 and measure time taken
batch_size = 10

println("Starting to parse each file and create dataframe")
@time begin
    data_frames = DataFrame[]  # Store parsed DataFrames
    for i in 1:batch_size:length(selected_files)
        batch_files = selected_files[i:min(i + batch_size - 1, end)]
        
        # Parse each file in the batch
        for file in batch_files
            data_df = parse_cfg_file(file)
            # Check if parsing succeeded
            if data_df === nothing
                println("Skipping $file due to parsing error: $file")
                continue
            end

            data_df = rename_columns(data_df)
            push!(data_frames, data_df)
        end

        # Insert parsed batch data
        if !isempty(data_frames)
            batch_data = vcat(data_frames...)
            insert_dataframe_to_postgresql(conn, batch_data, table_name)
            print(".")
            empty!(data_frames)
        end
    end
end

# Close the connection
close(conn)
