from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from scipy.interpolate import make_interp_spline
from scipy.ndimage.filters import gaussian_filter1d

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 13
timestep_1 = 2 
tau = 1000 * timestep_1
delQ = 0.05 #eV


list_val = [
#region start
[ -41.8372 ,205.55 ],
[ -41.7371 ,233.204 ],
[ -41.6371 ,269.283 ],
[ -41.5371 ,300.264 ],
[ -41.437 ,338.188 ],
[ -41.3368 ,351.69 ],
[ -41.2358 ,359.911 ],
[ -41.1372 ,399.771 ],
[ -41.0356 ,413.256 ],
[ -40.9318 ,430.081 ],
[ -40.8288 ,395.459 ],
[ -40.7309 ,448.048 ],
[ -40.6332 ,518.031 ],
[ -40.533 ,529.319 ],
[ -40.431 ,591.617 ],
[ -40.3271 ,562.708 ],
[ -40.2239 ,564.759 ],
[ -40.1241 ,576.268 ],
[ -40.0225 ,603.968 ],
[ -39.9246 ,648.777 ],
[ -39.8202 ,614.398 ],
[ -39.723 ,668.86 ],
[ -39.6138 ,637.458 ],
[ -39.5146 ,645.288 ],
[ -39.4208 ,769.693 ],
[ -39.3149 ,712.827 ],
[ -39.2145 ,762.571 ],
[ -39.1123 ,692.985 ],
[ -39.0155 ,754.476 ],
[ -38.9152 ,866.674 ],
[ -38.8134 ,835.377 ],
[ -38.7078 ,914.012 ],
[ -38.6132 ,926.156 ],
[ -38.5191 ,1053.19 ],
[ -38.4109 ,842.22 ],
[ -38.3055 ,882.338 ],
[ -38.2139 ,925.13 ],
[ -38.1151 ,1097.02 ],
[ -38.0114 ,1000.53 ],

#endRegion
]

curve_list_1 = [
#region start
[ -41.8372 ,205.55 ],
[ -41.7371 ,233.204 ],
[ -41.6371 ,269.283 ],
[ -41.5371 ,300.264 ],
[ -41.437 ,338.188 ],
[ -41.3368 ,351.69 ],
[ -41.2358 ,359.911 ],
[ -41.1372 ,399.771 ],
[ -41.0356 ,413.256 ],
[ -40.9318 ,430.081 ],
[ -40.8288 ,395.459 ],
[ -40.7309 ,448.048 ],
[ -40.6332 ,518.031 ],
[ -40.533 ,529.319 ],
[ -40.431 ,591.617 ],

#endRegion
]

curve_list_2 = [
#region start

[ -39.5146 ,645.288 ],
[ -39.4208 ,769.693 ],
[ -39.3149 ,712.827 ],
[ -39.2145 ,762.571 ],
[ -39.1123 ,692.985 ],
[ -39.0155 ,754.476 ],
[ -38.9152 ,866.674 ],
[ -38.8134 ,835.377 ],
[ -38.7078 ,914.012 ],
[ -38.6132 ,926.156 ],
[ -38.5191 ,1053.19 ],
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

ax.text(0.25,0.75,f'Heat Capacity - {round(a/atoms_num,5)} eV/K \nMelting point - {melting} K \nLatent heat - {round((40.431 - 39.9246)/atoms_num,2)} eV' ,ha='center', va='center',fontsize=10,transform=ax.transAxes,  bbox=dict(facecolor='red', alpha=0.5) )
ax.legend(["Simulated","Curvefit"])
path = os.path.join(f"plot_code/cluster_sizes")
os.makedirs(path, exist_ok=True)
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
# fig.savefig(save_path, bbox_inches='tight')
plt.show()

# plt.text(x_axis[0],-0.079,f'Heat Capacity - {round(a,5)/atoms_num} eV/K \nMelting point - {melting} K \nLatent heat - {round((491.461- 485.395),2)} eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )

