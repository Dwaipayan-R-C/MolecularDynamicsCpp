from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from scipy.interpolate import make_interp_spline
from scipy.ndimage.filters import gaussian_filter1d

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 2869
timestep_1 = 2 
tau = 1000 * timestep_1
delQ = 50 #eV

list_val = [
#region start
[ -10447.5 ,189.971 ],
[ -10395 ,251.517 ],
[ -10342.2 ,288.424 ],
[ -10289.9 ,346.138 ],
[ -10240.8 ,409.142 ],
[ -10192.7 ,504.735 ],
[ -10143.2 ,588.783 ],
[ -10092 ,642.166 ],
[ -10041 ,673.529 ],
[ -9990.29 ,717.714 ],
[ -9939.71 ,765.122 ],
[ -9889.28 ,804.066 ],
[ -9839.05 ,846.273 ],
[ -9788.93 ,893.881 ],
[ -9738.56 ,930.59 ],
[ -9688.12 ,945.706 ],
[ -9637.7 ,951.686 ],
[ -9587.14 ,952.516 ],
[ -9536.73 ,980.791 ],
[ -9486.11 ,1003.17 ],
[ -9435.93 ,1062.83 ],
[ -9385.52 ,1106.03 ],
[ -9335.27 ,1150.7 ],
[ -9285.04 ,1208.83 ],
#endRegion
]


curve_list_1 = [
#region start
[ -10240.8 ,409.142 ],
[ -10192.7 ,504.735 ],
[ -10143.2 ,588.783 ],
[ -10092 ,642.166 ],
[ -10041 ,673.529 ],
[ -9990.29 ,717.714 ],
[ -9939.71 ,765.122 ],
[ -9889.28 ,804.066 ],
[ -9839.05 ,846.273 ],
[ -9788.93 ,893.881 ],


#endRegion
]


curve_list_2 = [
#region start
[ -9435.93 ,1062.83 ],
[ -9385.52 ,1106.03 ],
[ -9335.27 ,1150.7 ],
[ -9285.04 ,1208.83 ],
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
plt.plot(x_axis,y_axis, color = 'brown')
plt.scatter(x_axis,y_axis,  color = 'black')
# plt.plot(curve_fit_x, a*curve_fit_x+b, color='blue', linewidth=1.5)
# plt.plot(curve_fit_x_high, a_high*curve_fit_x_high+b_high, color='blue', linewidth=1.5)
plt.grid()
#endregion

plt.text(x_axis[0],-9600,f'Heat Capacity - {round(a,5)} eV/K \nMelting point - {melting} K \nLatent heat - {round((9738.56 - 9587.14),2)} eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )
plt.legend(["Simulated","Curvefit"])
path = os.path.join(f"plot_code/cluster_sizes")
os.makedirs(path, exist_ok=True)
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
plt.savefig(save_path, bbox_inches='tight')
plt.show()
