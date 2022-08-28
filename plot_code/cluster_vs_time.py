from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from scipy.interpolate import make_interp_spline
from scipy.ndimage.filters import gaussian_filter1d

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 923
timestep_1 = 2 
tau = 1000 * timestep_1
delQ = 150 #eV


list_val = [
#region start
[147,2.537],
[309,6.448],
[923,24.798],
[2869,84.306],
[3871,114.748],
[5083,154.980],
[8217,257.185],
[10179,327.900],

#endRegion
]
# plt.xscale('log')
# plt.yscale('log')
#region plot start
x_axis = Extract(list_val,0)
y_axis = Extract(list_val,1)
# plt.plot(xnew, y_smooth, color='red')
plt.plot(x_axis, y_axis,'r-o')
plt.xlabel("Cluster size")
plt.ylabel("Time (seconds)")
plt.title(f"Time vs Cluster size (10000 iterations)")

# plt.plot(x_axis,y_axis, color='crimson')
# plt.plot(x_axis,f2)
plt.grid(True,which="both")
#endregion

# Save fig
path = os.path.join(f"plot_code/milestone_plots")
save_path = os.path.join(path,f"Cluster_vs_time.png")
os.makedirs(path, exist_ok=True)
# plt.text(x_axis[0],-3100,f'Heat Capacity - {heat_cap} eV/K \nMelting point - 750 K \nLatent heat - 127.5 eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )
# plt.xticks(np.arange(min(x_axis), max(x_axis)+1, 200))
plt.savefig(save_path, bbox_inches='tight')
plt.show()
