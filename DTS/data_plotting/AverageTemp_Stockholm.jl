using HTTP
using JSON
using Dates

# Define the base URL and parameters
base_url = "https://archive-api.open-meteo.com/v1/archive"

# Stockholm coordinates
latitude = 59.3293
longitude = 18.0686

# Date range
start_date = "2020-01-01"
end_date = "2024-12-31"

# Choose the variable for average temperature
variables = "temperature_2m_mean"

# Construct the query
query = "?latitude=$(latitude)&longitude=$(longitude)&start_date=$(start_date)&end_date=$(end_date)&daily=$(variables)&timezone=auto"

# Make the API request
response = HTTP.get(base_url * query)

# Parse the JSON response
data = JSON.parse(String(response.body))

# Extract dates and average temperatures
dates = data["daily"]["time"]
average_temp = data["daily"]["temperature_2m_mean"]

# Combine dates and average temperatures into a dictionary
temperature_data = Dict("dates" => dates, "average_temperature" => average_temp)
temperature_dict = Dict(temperature_data["dates"][i] => temperature_data["average_temperature"][i] for i in 1:length(temperature_data["dates"]))

# Optionally save the data to a file
open("stockholm_avg_temperature.json", "w") do file
    write(file, JSON.json(temperature_dict))
end

AirT = JSON.parsefile("C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\DTS\\stockholm_avg_temperature.json")
lookup_date = "2021-10-15"  # Change to the date you want to check
get(AirT, lookup_date, missing)