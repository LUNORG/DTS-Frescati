using DataFrames, CSV
using Luxor
using Statistics
using Dierckx
using Colors
using Random

df_BH1 = DataFrame(CSV.read("C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\DTS\\monthly_Tave_BH1.csv", DataFrame))
df_BH10 = DataFrame(CSV.read("C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\DTS\\monthly_Tave_BH1.csv", DataFrame))
weekly_BH1 = DataFrame(CSV.read("C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\DTS\\weekly_Tave_BH1.csv", DataFrame))
weekly_BH10 = DataFrame(CSV.read("C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\DTS\\weekly_Tave_BH10.csv", DataFrame))

function nan_filling(df::DataFrame)
    x = vcat(Vector(df[1,2:end]), Vector(df[2,2:end]), Vector(df[3,2:end]), Vector(df[4,2:end]), Vector(df[5,2:end]))
    known_idx = findall(!isnan, x)
    missing_idx = findall(isnan, x)
    t = collect(1:length(x))
    spline = Spline1D(t[known_idx], x[known_idx], k=3)  # k=3 means cubic spline
    x[missing_idx] .= spline(t[missing_idx])
    # X_matrix = Matrix(reshape(x[1:(5 * 12)], 12, 5)')
    return x
end
weekly_anomalies = nan_filling(weekly_BH10)[1:end-10]

month_T = nan_filling(df_BH10)
months = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]

function weekly_draw_frame(k)
    width, height = 900, 900
    tenRadius = 112
    zeroRadius = 45
    monthRadius = 180
    years = repeat(2020:2024, inner=48) 
        
    Drawing(width, height, "output.png")  # Create a drawing surface
    background("black")  # Set background to black
    origin()  # Move to center of the canvas

    setcolor("yellow")
    setline(1)

    # First Circle
    circle(Point(0, 0), zeroRadius*2, :stroke)  # Diameter 100 → Radius 50
    setfont("Helvetica", 24)
    Luxor.settext("0°", Point(zeroRadius*2+14, 0); halign="center", valign="center")

    # Second Circle
    setcolor("green")
    circle(Point(0, 0), tenRadius*2, :stroke)
    Luxor.settext("10°", Point(tenRadius*2+17, 0); halign="center", valign="center")

    # Third Circle
    setcolor("yellow")
    circle(Point(0, 0), monthRadius*2, :stroke) 

    setcolor("white")
    setfont("Helvetica", 45)
    Luxor.settext("$(years[k])", Point(0, 0); halign="center", valign="center")

    setcolor("yellow")
    setfont("Helvetica", 27)
    for i in 1:length(months)
        angle = ((i - 1) / length(months)) * 2π - π/2
        x = (monthRadius*2+26) * cos(angle)
        y = (monthRadius*2+26) * sin(angle)

        Luxor.settext(months[i], Point(x, y); halign="center", valign="center")
    end

    # Extract the year and anomaly data from the DataFrame
    newpath()
    setline(2)  # Set stroke weight
    sethue("white")  # Ensure stroke is white

    # Loop through the months and plot the anomalies
    for i in 1:k
        anomaly = weekly_anomalies[i]  # Get anomaly for the current month
        cycle_index = (i - 1) % 48 + 1  # This resets every 12 steps
        angle = (cycle_index - 1) / 48 * 2π - π/2  # Calculate the angle
        # angle = (i - 1) / length(months) * 2π - π/3  # Calculate angle for each month
        r = rescale(anomaly, 0, 10, zeroRadius*2, tenRadius*2)  # Map anomaly to radius
        x = r * cos(angle)
        y = r * sin(angle)

        line(x,y)
    end
        
    strokepath()  # Draw the shape outline
    finish()
    preview()
end

weekly_draw_frame(230)

function monthly_draw_frame(k)
    width, height = 600, 600
    tenRadius = 50
    zeroRadius = 20
    monthRadius = 80
    years = repeat(2020:2024, inner=12) 
        Drawing(width, height, "output.png")  # Create a drawing surface
        background("black")  # Set background to black
        origin()  # Move to center of the canvas

        setcolor("white")
        setline(1)

        # First Circle
        circle(Point(0, 0), zeroRadius*2, :stroke)  # Diameter 100 → Radius 50
        setfont("Helvetica", 50)
        setcolor("white")
        text("0°", Point(zeroRadius*2+8, 0), halign=:center)

        # Second Circle
        circle(Point(0, 0), tenRadius*2, :stroke)
        text("10°", Point(tenRadius*2+8, 0), halign=:center)

        # Third Circle
        circle(Point(0, 0), monthRadius*2, :stroke) 
        
        text("$(years[k])", Point(0, 0), halign=:center, valign=:center)

        for i in 1:length(months)
            angle = ((i - 1) / length(months)) * 2π - π/3
            x = (monthRadius*2+14) * cos(angle)
            y = (monthRadius*2+14) * sin(angle)

            text(months[i], Point(x, y), halign=:center, valign=:middle)
        end

        # Extract the year and anomaly data from the DataFrame
        newpath()
        setline(2)  # Set stroke weight
        sethue("white")  # Ensure stroke is white

        # Loop through the months and plot the anomalies
        for i in 1:k
            anomaly = month_T[i]  # Get anomaly for the current month
            cycle_index = (i - 1) % 12 + 1  # This resets every 12 steps
            angle = (cycle_index - 1) / 12 * 2π - π/3  # Calculate the angle
            # angle = (i - 1) / length(months) * 2π - π/3  # Calculate angle for each month
            r = rescale(anomaly, 0, 10, zeroRadius*2, tenRadius*2)  # Map anomaly to radius
            x = r * cos(angle)
            y = r * sin(angle)

            line(x,y)
        end
        
        strokepath()  # Draw the shape outline
        finish()
        preview()
end

using Plots
using Images
using FFMPEG

anim = Animation()

anim = @animate for k in 1:230
    weekly_draw_frame(k)
    img = Images.load("C:\\Users\\matil\\OneDrive\\Documenti\\UNIBO\\master_thesis_KTH\\DTS\\data_analysis\\output.png")
    # Plot the image (you can customize the plotting as needed)
    heatmap(Matrix(img), color=:auto, axis=false, size=(900, 900), dpi=300)
    # heatmap(Matrix(img), color=:auto, axis=false)  # Adjust the plot style as needed
end

gif(anim, "climatespiral_BH10.gif", fps=10)
