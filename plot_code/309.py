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
delQ = 3.5 #eV


list_val = [
#region start
[ -1104.91 ,88.0259 ],
[ -1101.56 ,134.421 ],
[ -1097.8 ,173.308 ],
[ -1094.44 ,222.692 ],
[ -1090.88 ,267.041 ],
[ -1087.36 ,314.403 ],
[ -1083.74 ,327.166 ],
[ -1080.36 ,393.702 ],
[ -1076.7 ,425.65 ],
[ -1073.15 ,453.989 ],
[ -1069.71 ,476.384 ],
[ -1066.14 ,508.725 ],
[ -1062.65 ,562.403 ],
[ -1059.12 ,601.441 ],
[ -1055.62 ,643.072 ],
[ -1052.11 ,659.571 ],
[ -1048.58 ,708.814 ],
[ -1045.03 ,735.77 ],
[ -1041.53 ,761.225 ],
[ -1037.97 ,751.491 ],
[ -1034.44 ,760.366 ],
[ -1030.92 ,750.749 ],
[ -1027.39 ,766.474 ],
[ -1023.87 ,790.7 ],
[ -1020.33 ,829.778 ],
[ -1016.82 ,848.255 ],
[ -1013.29 ,870.812 ],
[ -1009.77 ,920.161 ],
[ -1006.22 ,947.139 ],
[ -1002.73 ,951.922 ],
[ -999.169 ,998.564 ],
[ -995.703 ,1043.64 ],
[ -992.172 ,1077.96 ],
[ -988.664 ,1089.31 ],
[ -985.127 ,1151.82 ],
[ -981.605 ,1175.69 ],
[ -978.041 ,1197.64 ],
[ -974.523 ,1223.36 ],

#endRegion
]

curve_list_1 = [
#region start

[-1080.36 ,393.702 ],
[ -1076.7 ,425.65 ],
[ -1073.15 ,453.989 ],
[ -1069.71 ,476.384 ],
[ -1066.14 ,508.725 ],
[ -1062.65 ,562.403 ],
[ -1059.12 ,601.441 ],
[ -1055.62 ,643.072 ],
[ -1052.11 ,659.571 ],
[ -1048.58 ,708.814 ],
[ -1045.03 ,735.77 ],
[ -1041.53 ,761.225 ],
#endRegion
]

curve_list_2 = [
#region start
[ -1016.82 ,848.255 ],
[ -1013.29 ,870.812 ],
[ -1009.77 ,920.161 ],
[ -1006.22 ,947.139 ],
[ -1002.73 ,951.922 ],
[ -999.169 ,998.564 ],
[ -995.703 ,1043.64 ],
[ -992.172 ,1077.96 ],
[ -988.664 ,1089.31 ],
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
plt.plot(curve_fit_x, a*curve_fit_x+b, color='blue', linewidth=1.5)
plt.plot(curve_fit_x_high, a_high*curve_fit_x_high+b_high, color='blue', linewidth=1.5)
plt.grid()
#endregion

path = os.path.join(f"{atoms_num} atoms")
save_path = os.path.join(path,f"{tau}_{timestep_1}_{delQ}_EvsT.png")
os.makedirs(path, exist_ok=True)
plt.text(x_axis[0],-1010,f'Heat Capacity - {round(a,5)} eV/K \nMelting point - {melting} K \nLatent heat - {round((1041.53 - 1030.92),2)} eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )
plt.legend(["Simulated","Curvefit"])
plt.show()
# # Save fig
path = os.path.join(f"plot_code/cluster_sizes")
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
# plt.savefig(save_path, bbox_inches='tight')
# plt.show()
