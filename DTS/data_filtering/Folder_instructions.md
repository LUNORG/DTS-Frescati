This folder is dedicated to filtering the data.

The optical fiber are longer than the boreholes and there are discradable parts. In order to evaluate which section of the fiber is actually within the boreholes, the metadata values do not align with reality (they're probably just default values). 
The filtering method adopted includes 2 main criteria: the absolute T is set to be within a chosen range (-5°C-30°C), and the difference between two consecutive (in space) tempertaure values has to be lower than a chosen threshold (1.0 °C).

The longest segment where this criteria are satisfied for more than 90% of the time (I considered all the available files), is the right one.

Borehole 10 is 350m long, while Borehole 1 is 100m long. Knowing this, from documentation, I hold accountable the last point from the previous filtering method, while the beginning point would be obtained substracting the boreholes' length.

Additional filtering is applied using a Kalman filter that is based on the variance value obtained through the uncertanty evaluation procedure. In this case the purpose of filtering was reducing the noise, smoothening the data.