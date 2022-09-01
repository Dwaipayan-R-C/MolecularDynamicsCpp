from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from scipy.interpolate import make_interp_spline
from scipy.ndimage.filters import gaussian_filter1d

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 5083
timestep_1 = 2 
tau = 1000 * timestep_1
delQ = 120 #eV


list_val = [
#region start
[ -18549.7 ,195.464 ],
[ -18430.5 ,278.67 ],
[ -18305.7 ,362.892 ],
[ -18183.4 ,426.585 ],
[ -18061.6 ,509.926 ],
[ -17940.8 ,567.441 ],
[ -17819.7 ,629.458 ],
[ -17699.6 ,733.411 ],
[ -17580.4 ,828.008 ],
[ -17460 ,884.727 ],
[ -17338.6 ,913.853 ],
[ -17217.2 ,929.497 ],
[ -17096.5 ,958.793 ],
[ -16975.9 ,992.561 ],
[ -16855.1 ,1032.39 ],
[ -16734.6 ,1094.53 ],
[ -16613.4 ,1150.42 ],
[ -16492.8 ,1220.92 ],
[ -16372.2 ,1302.49 ],
[ -16251 ,1352.97 ],
[ -16130.1 ,1418.83 ],
[ -16009.3 ,1490.43 ],


#endRegion
]


curve_list_1 = [
#region start
[ -18549.7 ,195.464 ],
[ -18430.5 ,278.67 ],
[ -18305.7 ,362.892 ],
[ -18183.4 ,426.585 ],
[ -18061.6 ,509.926 ],
[ -17940.8 ,567.441 ],
[ -17819.7 ,629.458 ],
[ -17699.6 ,733.411 ],
[ -17580.4 ,828.008 ],
[ -17460 ,884.727 ],

#endRegion
]


curve_list_2 = [
#region start
[ -16975.9 ,992.561 ],
[ -16855.1 ,1032.39 ],
[ -16734.6 ,1094.53 ],
[ -16613.4 ,1150.42 ],
[ -16492.8 ,1220.92 ],
[ -16372.2 ,1302.49 ],
[ -16251 ,1352.97 ],
[ -16130.1 ,1418.83 ],
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
# ax.plot(curve_fit_x, a*curve_fit_x+b, color='blue', linewidth=1.5)
# ax.plot(curve_fit_x_high, a_high*curve_fit_x_high+b_high, color='blue', linewidth=1.5)
ax.grid()
#endregion

ax.text(0.3,0.75,f'Heat Capacity - {round(a/atoms_num,6)} eV/(K*atoms) \nMelting point - {melting} K \nLatent heat - {round((17460 - 16975.9 )/atoms_num,4)} eV/atom' ,ha='center', va='center',fontsize=10,transform=ax.transAxes,  bbox=dict(facecolor='red', alpha=0.5) )
ax.legend(["Simulated","Curvefit"])
path = os.path.join(f"plot_code/cluster_sizes")
os.makedirs(path, exist_ok=True)
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
fig.savefig(save_path, bbox_inches='tight')
plt.show()
# plt.text(x_axis[0],-16800,f'Heat Capacity - {round(a/atoms_num,5)} eV/K \nMelting point - {melting} K \nLatent heat - {round((17460 - 16975.9),2)} eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )

