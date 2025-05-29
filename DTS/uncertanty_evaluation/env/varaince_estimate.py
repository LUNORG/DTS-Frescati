# set the python interpreter as & C:/Users/matil/.conda/envs/dtscalibration_env/python.exe 
# on the terminal type:
# $env:PYTHONPATH="C:\Users\matil\.conda\envs\dtscalibration_env\Lib\site-packages"
import json
import os
import glob
from dtscalibration import read_sensornet_files
from dtscalibration.variance_stokes import variance_stokes_constant

# Function to process a batch of files and calculate variance
def process_batch(filepath, Lmin, Lmax):
    ds = read_sensornet_files(directory=filepath)
    ds = ds.sel(x=slice(Lmin, Lmax))
    sections = { "referenceTemperature": [slice(-38, -5)]}
    st_var, resid = variance_stokes_constant(
        ds.dts.st, sections, ds.dts.acquisitiontime_fw, reshape_residuals=True
    )
    ast_var, _ = variance_stokes_constant(
        ds.dts.ast, sections, ds.dts.acquisitiontime_fw, reshape_residuals=False
    )
    out = ds.dts.calibrate_single_ended(sections=sections, st_var=st_var, ast_var=ast_var)
    ds1 = out.isel(time=0)
    stast_var = ds1.var_fw_da.sel(comp_fw=["dT_dst", "dT_dast"]).sum(dim="comp_fw")
    return stast_var.values, ds1.x.values  # Return both variance and length values

filepath = 'E:\\partial channel 2'

# channel1_variances, channel1_lengths = process_batch(filepath, 40, 160) but it has to start from -38
# channel2_variances, channel2_lengths = process_batch(filepath, 0, 412) but it has to start from -38

channel2_variances, channel2_lengths = process_batch(filepath, -38, 160)

# variance_dict_channel1 = {length: variance for length, variance in zip(channel1_lengths, channel1_variances)}
variance_dict_channel2 = {length: variance for length, variance in zip(channel2_lengths, channel2_variances)}

# Save dictionary to a JSON file
with open("variance_dict_channel2.json", "w") as f:
    json.dump(variance_dict_channel2, f)