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
mytime = 21
vmax = 0.5#np.max(abs(data[mytime,:,:]*20))
bins = np.linspace(-vmax, vmax, 10, endpoint=True)
k = data[mytime]


super_threshold_indices = np.abs(k) < 1e-6
k[super_threshold_indices] = 0
mydata = np.digitize(data[mytime]*20,bins)

cax = plt.contourf(k*20,levels=bins, cmap="seismic", vmin=-vmax, vmax=vmax, extend="both", norm="logit")
#plt.contour(np.array(range(512))/8, np.array(range(512))/8, data[mytime]*20,levels=bins, colors="gray", linewidths=0.5, linestyles="solid", negative_linestyles="solid")
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
print(ds["time"][mytime])
print(super_threshold_indices)

