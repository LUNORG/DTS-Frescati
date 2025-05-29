This folder is dedicated to the creation of various analysis.

*distance_distribution.jl* is generated with pluto.jl and recreates the 3D structure of the borehole field. The code uses the deviation measurements available, regarding the inclination of each borehole, and the exact position in space was recreated using design information on https://www.sgu.se/en/products/maps/map-viewer/groundwater-map-viewers/wells/ 

*sin_fitting.jl* fits for the 3 datasets BH1, BH10 and the outdoor temperature a sinusoidal function in the form $A+Bsin(2\pi/T+phi) to the raw data. Observation regarding phase shift, amplitude damping can be made.

*undisturbedT.jl* this code uses the temperature data provided by the machine used for the deviation measurements, therefore relative to 2 days after the boreholes drilling and before they were utilized and connected. This undisturbed data can be compared to the data provided by DTS.

*Fourier.jl* based on FFTW package a discrete fast fourier analysis was made for the 3 datasets, evaluating frequency, corresponding amplitude and phase. The dominant periodicity was the same and equal to 350days.

*monthlyAve.jl* code for computing montlhy, weekly and daily temperature averages.

*pointsource.jl* code to model a point heat source in an infinite homogeneous solid. The load is the sum of multiple step functions that model a sinusoidal load. The goal is to evaluate the temperature at various distances from it.