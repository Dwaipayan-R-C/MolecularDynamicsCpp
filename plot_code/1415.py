from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from scipy.interpolate import make_interp_spline
from scipy.ndimage.filters import gaussian_filter1d

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 1415
timestep_1 = 2 
tau = 1000 * timestep_1
delQ = 15 #eV

list_val = [
#region start
[ -5131.94 ,143.972 ],
[ -5116.22 ,179.148 ],
[ -5100.61 ,223.602 ],
[ -5086.02 ,280.137 ],
[ -5069.99 ,303.705 ],
[ -5055.53 ,360.061 ],
[ -5039.78 ,378.933 ],
[ -5025.54 ,431.531 ],
[ -5009.83 ,450.993 ],
[ -4995.09 ,495.263 ],
[ -4979.57 ,518.102 ],
[ -4964.68 ,569.643 ],
[ -4949.64 ,594.299 ],
[ -4934.51 ,633.741 ],
[ -4919.4 ,663.636 ],
[ -4904.3 ,698.544 ],
[ -4889.22 ,736.355 ],
[ -4874.04 ,754.352 ],
[ -4858.96 ,787.058 ],
[ -4843.84 ,804.816 ],
[ -4828.74 ,846.884 ],
[ -4813.59 ,857.844 ],
[ -4798.48 ,879.983 ],
[ -4783.38 ,884.176 ],
[ -4768.26 ,898.552 ],
[ -4753.12 ,888.098 ],
[ -4738.08 ,896.719 ],
[ -4722.84 ,888.446 ],
[ -4707.84 ,930.459 ],
[ -4692.75 ,958.321 ],
[ -4677.57 ,992.929 ],
[ -4662.46 ,1008.14 ],
[ -4647.35 ,1042.44 ],
[ -4632.31 ,1072.27 ],
[ -4617.23 ,1092.14 ],
[ -4602.01 ,1105.94 ],
[ -4587.04 ,1166.34 ],
[ -4571.75 ,1165.66 ],
[ -4556.8 ,1219.97 ],


#endRegion
]


curve_list_1 = [
#region start
[ -5009.83 ,450.993 ],
[ -4995.09 ,495.263 ],
[ -4979.57 ,518.102 ],
[ -4964.68 ,569.643 ],
[ -4949.64 ,594.299 ],
[ -4934.51 ,633.741 ],
[ -4919.4 ,663.636 ],
[ -4904.3 ,698.544 ],
[ -4889.22 ,736.355 ],
[ -4874.04 ,754.352 ],
[ -4858.96 ,787.058 ],
[ -4843.84 ,804.816 ],
[ -4828.74 ,846.884 ],


#endRegion
]


curve_list_2 = [
#region start
[ -4692.75 ,958.321 ],
[ -4677.57 ,992.929 ],
[ -4662.46 ,1008.14 ],
[ -4647.35 ,1042.44 ],
[ -4632.31 ,1072.27 ],
[ -4617.23 ,1092.14 ],
[ -4602.01 ,1105.94 ],
[ -4587.04 ,1166.34 ],
[ -4571.75 ,1165.66 ],
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

plt.text(x_axis[0],-4750,f'Heat Capacity - {round(a,5)} eV/K \nMelting point - {melting} K \nLatent heat - {round((4828.74 - 4722.84),2)} eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )
plt.legend(["Simulated","Curvefit"])
path = os.path.join(f"plot_code/cluster_sizes")
os.makedirs(path, exist_ok=True)
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
plt.savefig(save_path, bbox_inches='tight')
plt.show()
