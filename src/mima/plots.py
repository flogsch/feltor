import matplotlib.pyplot as plt
import netCDF4 as nc 
import json
import numpy as np

with open("input/default.json") as f:
    d = json.load(f)
    #print(d)
Nx = d["output"]["Nx"]
Ny = d["output"]["Ny"]
#print(Nx)
fn = "mima.nc"
ds = nc.Dataset(fn)
name = "vort"
# print(ds)
# print(ds["phi"][:,:,1])
data = ds[name][50]
mydata = np.flip(data, axis=1)
vmax = np.max(abs(data))
plt.pcolormesh(data, cmap="bwr", vmin=-vmax, vmax=vmax)
plt.colorbar()
plt.title(name)
plt.show() 