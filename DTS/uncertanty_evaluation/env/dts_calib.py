# set the python interpreter as & C:/Users/matil/.conda/envs/dtscalibration_env/python.exe 
# on the terminal type:
# $env:PYTHONPATH="C:\Users\matil\.conda\envs\dtscalibration_env\Lib\site-packages"

import os
import glob

from dtscalibration import read_sensornet_files
from dtscalibration.variance_stokes import variance_stokes_constant
import matplotlib.pyplot as plt

# Set the directory containing your .ddf files
filepath = '....'
filepathlist = sorted(glob.glob(os.path.join(filepath, "*.ddf")))
filenamelist = [os.path.basename(path) for path in filepathlist]

# for fn in filenamelist: print(fn)

ds = read_sensornet_files(directory=filepath)
# print(ds)

plt.ion()

ds = ds.sel(x=slice(-38, 412))  # dismiss parts of the fiber that are not interesting
                                # (40,160) for channel 1 and (0, 412) for channel 2 but it must include the refernce section

sections = {
    "referenceTemperature": [slice(-38, -5)]  # Use reference temperature across the region
}

# sections = {
#   "probe1Temperature": [slice(10, 25.5)],  # warm bath
#    "probe2Temperature": [slice(-5.5, 5.5)],  # cold bath
#}

st_var, resid = variance_stokes_constant(
    ds.dts.st, sections, ds.dts.acquisitiontime_fw, reshape_residuals=True
)
ast_var, _ = variance_stokes_constant(
    ds.dts.ast, sections, ds.dts.acquisitiontime_fw, reshape_residuals=False
)

out = ds.dts.calibrate_single_ended(sections=sections, st_var=st_var, ast_var=ast_var)

ds1 = out.isel(time=0)

# Uncertainty from the noise in (anti-) stokes measurements
stast_var = ds1.var_fw_da.sel(comp_fw=["dT_dst", "dT_dast"]).sum(dim="comp_fw")

# Parameter uncertainty
par_var = ds1.tmpf_var - stast_var

# Plot
plt.figure(figsize=(12, 4))
plt.fill_between(ds1.x, stast_var, label="Noise in (anti-) stokes measurements")
plt.fill_between(ds1.x, y1=ds1.tmpf_var, y2=stast_var, label="Parameter estimation")
plt.suptitle("Variance of the estimated temperature")
plt.ylabel("$\\sigma^2$ ($^\\circ$C$^2$)")
plt.legend(fontsize="small")

# input("Press Enter to continue...")

# The effects of the parameter uncertainty can be further inspected
# Note that the parameter uncertainty is not constant over the fiber and certain covariations can reduce to temperature uncertainty
ds1.var_fw_da.plot(hue="comp_fw", figsize=(12, 4))

out2 = ds.dts.monte_carlo_single_ended(
    result=out,
    st_var=st_var,
    ast_var=ast_var,
    conf_ints=[2.5, 97.5],
    mc_sample_size=500,
)

ds1 = out.isel(time=0)

(out2.isel(time=0).tmpf_mc_var ** 0.5).plot(
    figsize=(12, 4), label="Monte Carlo approx."
)
(out.isel(time=0).tmpf_var ** 0.5).plot(label="Linear error approx.")
plt.ylabel("$\\sigma$ ($^\\circ$C)")
plt.legend(fontsize="small")

out.isel(time=0).tmpf.plot(linewidth=0.7, figsize=(12, 4))
out2.isel(time=0).tmpf_mc.sel(CI=2.5).plot(linewidth=0.7, label="CI: 2.5%")
out2.isel(time=0).tmpf_mc.sel(CI=97.5).plot(linewidth=0.7, label="CI: 97.5%")
plt.legend(fontsize="small")

input("Press enter to continue...")