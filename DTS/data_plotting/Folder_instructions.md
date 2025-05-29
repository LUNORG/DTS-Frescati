This folder is dedicated to the creation of various plots.

Every plot starts with spline objects, used for smoothening and interpolation purposes. One has to pick the rigth value of the smoothness parameter s. In order to do so an iterative method is the best approach: balancing the error introduced with high s and the noise effect.

*plotting_process.jl* returns contour, surface plots and heatmaps. It also generates contour plots of the difference between the temperature at each position and the average value over the length.  
*AveTemp_stockholm.jl* is dedicated to retrieving the average dialy temperature data from API open meteo website and storing it into a dictionary.
*watertable.jl* used for plotting a line at the supposed ground water level (-43m).
*Tave_plot.jl* used for plotting the average temperature over length in time, for the two boreholes and the outdoor temperature in stockholm as well, to enable a visual comparison.
*ClimateSpiral.jl* uses depth and weekly averaged temperature values, interpolates them using splines and plots them in an animation as a climat spiral.

LN: Our plots would include the .txt files right?