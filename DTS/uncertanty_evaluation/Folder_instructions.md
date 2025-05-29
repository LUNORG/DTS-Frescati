This folder is dedicated to the uncertanty evaluation

For this purpose dtscalibration python library was used. It provides functions to evaluate the temperature variance associated with each length. It requests ddf files as inputs with the same time step. I had to delete few files where the time step was not exaclty 300.05 seconds but 0.06 more or less than that. The higher the number of files provided to the function the better the more accurate the evaluation. 

To know more about the mentioned python package visit: "https://python-dts-calibration.readthedocs.io/en/stable/"