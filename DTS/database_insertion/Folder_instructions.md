This folder is dedicated to the insertion of DTS data in a PostgreSQL database. 

The data are in the form of ddf files, or cfg files. 
Each ddf file has a first section with metadata, followed by actual temperature data, associated with length, stokes and anti-stokes corresponding values. 

Each line would be parsed and the 2 sections are separated into 2 different dataframes. The metadata results in one single row per ddf file, while the data contains as many rows as the selected number of points where the fiber sensor takes measurements.

Each dataframe will be inserted into 2 different tables in the database, one DTSold_metadata and 1 for the data: DTSold_data.

They have different column haeders but can be linked to one another (data and metadata obtained from the same ddf file) thanks to the datetime and channel column. They identify the exact date and time of the measurement and creation of file, while channel is either 1 or 2: channel 1 referes to the fiber in Borehole 1, channel 2 refers to Borehole 10.

A preliminary selection was carried out because a large portion of the initial and final temperature values are completely unreliable (such as 700°C) and do not correspond to an actual part of the optical fiber. Therefore, only the central lengths were considered. For the identification of the limits, which are nevertheless very broad and not stringent, see the files: *ddf_max_length.jl* and *ddf_min_length.jl*

LN: so cfg isn't used at all? 
Our flow rate and temp in and out data are in .txt-files.