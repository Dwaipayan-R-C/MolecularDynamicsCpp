from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from scipy.interpolate import make_interp_spline
from scipy.ndimage.filters import gaussian_filter1d

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 55
timestep_1 = 2 
tau = 1000 * timestep_1
delQ = 0.3 #eV


list_val = [
#region start
[ -190.272 ,53.3126 ],
[ -190.104 ,72.3329 ],
[ -189.97 ,95.5578 ],
[ -189.794 ,115.662 ],
[ -189.63 ,132.496 ],
[ -189.502 ,160.079 ],
[ -189.339 ,178.599 ],
[ -189.155 ,195.43 ],
[ -188.916 ,205.976 ],
[ -188.779 ,232.171 ],
[ -188.565 ,242.659 ],
[ -188.584 ,286.597 ],
[ -188.149 ,271.261 ],
[ -188.128 ,310.781 ],
[ -187.993 ,333.953 ],
[ -187.849 ,356.192 ],
[ -187.507 ,352.052 ],
[ -187.467 ,387.96 ],
[ -187.299 ,406.502 ],
[ -187.187 ,434.146 ],
[ -186.974 ,446.153 ],
[ -186.832 ,468.924 ],
[ -186.515 ,468.416 ],
[ -186.289 ,479.321 ],
[ -186.045 ,488.579 ],
[ -186.102 ,538.722 ],
[ -185.711 ,525.674 ],
[ -185.457 ,532.923 ],
[ -185.364 ,561.483 ],
[ -185.077 ,564.172 ],
[ -184.747 ,561.576 ],
[ -184.482 ,566.069 ],
[ -184.254 ,576.513 ],
[ -183.953 ,577.582 ],
[ -183.979 ,623.085 ],
[ -183.695 ,626.059 ],
[ -183.762 ,676.343 ],
[ -183.568 ,693.354 ],
[ -183.274 ,693.047 ],
[ -182.879 ,680.463 ],
[ -182.962 ,734.251 ],
[ -182.629 ,729.519 ],
[ -182.407 ,741.504 ],
[ -181.999 ,727.284 ],
[ -182.046 ,775.51 ],
[ -181.859 ,790.492 ],


#endRegion
]

#region plot start
x_axis = Extract(list_val,1)
y_axis = Extract(list_val,0)
param = np.linspace(0, 1, len(x_axis))
spl = make_interp_spline(param, np.c_[x_axis,y_axis], k=2) #(1)
xnew, y_smooth = spl(np.linspace(0, 1, len(x_axis) * 100)).T #(2)
plt.plot(xnew, y_smooth, color='red')
plt.scatter(x_axis, y_axis, color="black")
plt.xlabel("Temperature (K)")
plt.ylabel("Total Energy (eV)")
plt.suptitle(f"Total Energy vs Temperature ( {atoms_num} atoms )")
plt.title(f"Tau = {tau} fs, timestep = {timestep_1} fs, ΔQ = {delQ} eV")
# plt.plot(x_axis,y_axis, color='crimson')
# plt.plot(x_axis,f2)

plt.grid()
#endregion

# Save fig

path = os.path.join(f"plot_code/cluster_sizes")
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
os.makedirs(path, exist_ok=True)
# plt.text(x_axis[0],-3100,f'Heat Capacity - {heat_cap} eV/K \nMelting point - 750 K \nLatent heat - 127.5 eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )
# plt.xticks(np.arange(min(x_axis), max(x_axis)+1, 200))
# plt.savefig(save_path, bbox_inches='tight')
plt.show()
