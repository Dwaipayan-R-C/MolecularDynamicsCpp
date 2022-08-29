from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from scipy.interpolate import make_interp_spline
from scipy.ndimage.filters import gaussian_filter1d

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 309
timestep_1 = 2 
tau = 1000 * timestep_1
delQ = 4 #eV


list_val = [
#region start
[ -1105.04 ,86.5576 ],
[ -1101.18 ,139.203 ],
[ -1096.92 ,184.091 ],
[ -1093.08 ,240.282 ],
[ -1088.97 ,281.633 ],
[ -1085.01 ,349.671 ],
[ -1080.94 ,363.817 ],
[ -1077.03 ,431.303 ],
[ -1072.93 ,463.795 ],
[ -1068.9 ,509.972 ],
[ -1064.87 ,571.685 ],
[ -1060.79 ,593.645 ],
[ -1056.83 ,645.325 ],
[ -1052.8 ,666.445 ],
[ -1048.76 ,702.247 ],
[ -1044.75 ,756.288 ],
[ -1040.69 ,747.027 ],
[ -1036.65 ,776.351 ],
[ -1032.62 ,791.788 ],
[ -1028.58 ,779.601 ],
[ -1024.56 ,831.366 ],
[ -1020.52 ,843.438 ],
[ -1016.48 ,875.407 ],
[ -1012.44 ,892.341 ],
[ -1008.42 ,931.419 ],
[ -1004.41 ,969.539 ],
[ -1000.36 ,1006.14 ],
[ -996.339 ,1020.02 ],
[ -992.325 ,1092.05 ],
[ -988.262 ,1110.63 ],


#endRegion
]

curve_list_1 = [
#region start
[ -1105.04 ,86.5576 ],
[ -1101.18 ,139.203 ],
[ -1096.92 ,184.091 ],
[ -1093.08 ,240.282 ],
[ -1088.97 ,281.633 ],
[ -1085.01 ,349.671 ],
[ -1080.94 ,363.817 ],
[ -1077.03 ,431.303 ],
[ -1072.93 ,463.795 ],
[ -1068.9 ,509.972 ],
[ -1064.87 ,571.685 ],
[ -1060.79 ,593.645 ],
[ -1056.83 ,645.325 ],
[ -1052.8 ,666.445 ],
[ -1048.76 ,702.247 ],
#endRegion
]

curve_list_2 = [
#region start
[ -1020.52 ,843.438 ],
[ -1016.48 ,875.407 ],
[ -1012.44 ,892.341 ],
[ -1008.42 ,931.419 ],
[ -1004.41 ,969.539 ],
[ -1000.36 ,1006.14 ],
[ -996.339 ,1020.02 ],
[ -992.325 ,1092.05 ],
[ -988.262 ,1110.63 ],



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

plt.xlabel("Temperature (K)")
plt.ylabel("Total Energy (eV)")
plt.suptitle(f"Total Energy vs Temperature ( {atoms_num} atoms )")
plt.title(f"Tau = {tau} fs, timestep = {timestep_1} fs, ΔQ = {delQ} eV")
# plt.legend(f"{timestep_1} fs", loc = 2)
# plt.plot(smoothed_mode(floor_range), floor_range,ls='solid', color='crimson')

plt.plot(x_axis,y_axis, color = 'brown')
plt.scatter(x_axis,y_axis,  color = 'black')
# plt.plot(curve_fit_x, a*curve_fit_x+b, color='blue', linewidth=1.5)
# plt.plot(curve_fit_x_high, a_high*curve_fit_x_high+b_high, color='blue', linewidth=1.5)
plt.grid()
#endregion

plt.text(x_axis[0],-1030,f'Heat Capacity - {round(a,5)} eV/K \nMelting point - {melting} K \nLatent heat - {round((1044.75- 1028.58),2)} eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )
plt.legend(["Simulated","Curvefit"])
# plt.show()
# # Save fig
path = os.path.join(f"plot_code/cluster_sizes")
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
plt.savefig(save_path, bbox_inches='tight')
# plt.show()
