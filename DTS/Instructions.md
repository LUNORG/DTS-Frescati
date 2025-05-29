This folder contains the code related to the processing of data coming from the DTS system.  

- **database_insertion**: Describes the process of copying all data to a local PC and subsequently inserting it into tables within a PostgreSQL database.  
- **data_filtering**: Focuses on the selection criteria for identifying the section of the optical fiber that is actually positioned inside the borehole, as well as the application of the Kalman filter.  
- **data_plotting**: Includes all functions for refining and smoothing the data, up to spline fitting, which enables visualization through various types of graphs.  
- **data_analysis**: Covers a range of different analyses performed on the data.  
- **uncertanty_evaluation**: Contains python files to compute the temperature variance of the data, using a specific python library.  

Each of these folders contains a **.md file** that provides a more detailed description of its contents.