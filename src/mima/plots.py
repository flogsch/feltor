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
name = "phi"
# print(ds)


tend = len(ds["time"])-1
data = ds[name]

#mydata = np.flip(data, axis=1)
mytime = 0
vmax = 5#np.max(abs(data[mytime,:,:]*20))
bins = np.linspace(-vmax, vmax, 10, endpoint=True)
mydata = np.digitize(data[mytime]*20,bins)

cax = plt.contourf(np.array(range(256))/4, np.array(range(256))/4, data[mytime]*20,levels=bins, cmap="seismic", vmin=-vmax, vmax=vmax, extend="both")
plt.contour(np.array(range(256))/4, np.array(range(256))/4, data[mytime]*20,levels=bins, colors="gray", linewidths=0.5, linestyles="solid", negative_linestyles="solid")
plt.colorbar(cax)

'''
fig,ax = plt.subplots()
cax = plt.pcolormesh(data[mytime]*20, cmap='seismic',vmin=-vmax, vmax=vmax)#, vmin=-vmax, vmax=vmax)
fig.colorbar(cax)
def animate(i):
    cax.set_array(data[i])
    plt.title(i)
#anim = FuncAnimation(fig, animate, interval=0.001, frames=400)
#anim.save('517.gif')
'''
plt.show()
print(ds["time"][1])
