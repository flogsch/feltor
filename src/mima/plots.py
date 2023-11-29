import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
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


tend = len(ds["time"])-1
data = ds[name]
#mydata = np.flip(data, axis=1)
vmax = np.max(abs(data[:,:,:]))

fig,ax = plt.subplots()
cax = plt.pcolormesh(data[0], cmap='bwr', vmin=-vmax/10, vmax=vmax/10)
fig.colorbar(cax)
def animate(i):
    cax.set_array(data[i])
    plt.title(i)
anim = FuncAnimation(fig, animate, interval=0.01, frames=10)
anim.save('517.gif')

plt.show()

