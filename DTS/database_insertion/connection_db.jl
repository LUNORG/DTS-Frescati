using LibPQ
using IniFile

# config.ini is the file containing the credentials needed for access
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

# Set conn to nothing initially
conn = nothing

# Attempt to establish a connection
try
    conn = LibPQ.Connection(conn_string)
    result = execute(conn, "SELECT version();")
    println("Connected to PostgreSQL successfully.")
    for row in result
        println(row)
    end
catch e
    println("Failed to connect to PostgreSQL: $e")
finally
    if conn !== nothing  # Check if conn was successfully created before closing
        close(conn)       # Always close the connection when done
    end
end
