from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from scipy.interpolate import make_interp_spline
from scipy.ndimage.filters import gaussian_filter1d

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 10179
timestep_1 = 2 
tau = 1000 * timestep_1
delQ = 190 #eV


list_val = [
#region start
[ -37212.3 ,180.794 ],
[ -37019.5 ,240.77 ],
[ -36839.9 ,349.801 ],
[ -36655.6 ,445.037 ],
[ -36458.9 ,437.338 ],
[ -36269.3 ,561.79 ],
[ -36078.5 ,613.83 ],
[ -35880.4 ,634.72 ],
[ -35690.3 ,692.367 ],
[ -35505.3 ,812.561 ],
[ -35310.9 ,832.748 ],
[ -35117.5 ,845.217 ],
[ -34928.6 ,925.933 ],
[ -34737.5 ,954.244 ],
[ -34545.2 ,959.574 ],
[ -34354 ,989.502 ],
[ -34163 ,1006.6 ],
[ -33971.2 ,1026.85 ],
[ -33779.9 ,1061.55 ],
[ -33588.5 ,1102.1 ],
[ -33397.4 ,1154.98 ],
[ -33206.4 ,1209.91 ],
[ -33014.9 ,1255.51 ],
[ -32823.6 ,1311.06 ],
[ -32631.8 ,1356.7 ],



#endRegion
]

curve_list_1 = [
#region start
[ -37212.3 ,180.794 ],
[ -37019.5 ,240.77 ],
[ -36839.9 ,349.801 ],
[ -36655.6 ,445.037 ],
[ -36458.9 ,437.338 ],
[ -36269.3 ,561.79 ],
[ -36078.5 ,613.83 ],
[ -35880.4 ,634.72 ],
[ -35690.3 ,692.367 ],
[ -35505.3 ,812.561 ],
[ -35310.9 ,832.748 ],
[ -35117.5 ,845.217 ],
[ -34928.6 ,925.933 ],


#endRegion
]

curve_list_2 = [
#region start
[ -33971.2 ,1026.85 ],
[ -33779.9 ,1061.55 ],
[ -33588.5 ,1102.1 ],
[ -33397.4 ,1154.98 ],
[ -33206.4 ,1209.91 ],
[ -33014.9 ,1255.51 ],
[ -32823.6 ,1311.06 ],
[ -32631.8 ,1356.7 ],

#endRegion
]

melting = np.array(Extract(curve_list_1,1))[-1]
#region plot start
x_axis = Extract(list_val,1)
y_axis = Extract(list_val,0)

curve_fit_x = np.array(Extract(curve_list_1,1))
curve_fit_y = np.array(Extract(curve_list_1,0))
curve_fit_x_high = np.array(Extract(curve_list_2,1))
curve_fit_y_high = np.array(Extract(curve_list_2,0))

a, b = np.polyfit(curve_fit_x, curve_fit_y, 1)
a_high, b_high = np.polyfit(curve_fit_x_high, curve_fit_y_high, 1)

fig,ax = plt.subplots()
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Total Energy (eV)")
fig.suptitle(f"Total Energy vs Temperature ( {atoms_num} atoms )")
ax.set_title(f"Tau = {tau} fs, timestep = {timestep_1} fs, ΔQ = {delQ} eV")
ax.plot(x_axis,y_axis, color = 'brown')
ax.scatter(x_axis,y_axis,  color = 'black')
ax.plot(curve_fit_x, a*curve_fit_x+b, color='blue', linewidth=1.5)
ax.plot(curve_fit_x_high, a_high*curve_fit_x_high+b_high, color='blue', linewidth=1.5)
ax.grid()
#endregion

ax.text(0.3,0.75,f'Heat Capacity - {round(a/atoms_num,6)} eV/(K*atoms) \nMelting point - {melting} K \nLatent heat - {round((34928.6-33971.2)/atoms_num,4)} eV/atom' ,ha='center', va='center',fontsize=10,transform=ax.transAxes,  bbox=dict(facecolor='red', alpha=0.5) )
ax.legend(["Simulated","Curvefit"])
path = os.path.join(f"plot_code/cluster_sizes")
os.makedirs(path, exist_ok=True)
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
# fig.savefig(save_path, bbox_inches='tight')
plt.show()
# plt.text(x_axis[0],-34000,f'Heat Capacity - {round(a/atoms_num,5)} eV/K \nMelting point - {melting} K \nLatent heat - {round((34728.1-33961.8),2)} eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )
