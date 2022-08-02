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
delQ = 1.5 #eV


list_val = [
#region start
[ -1108.42 ,88.0259 ],
[ -1107.93 ,109.111 ],
[ -1106.86 ,124.049 ],
[ -1106.32 ,146.785 ],
[ -1105.69 ,166.985 ],
[ -1104.53 ,178.909 ],
[ -1104.25 ,202.41 ],
[ -1102.62 ,206.679 ],
[ -1102.86 ,244.688 ],
[ -1101.33 ,250.178 ],
[ -1100.51 ,268.247 ],
[ -1100.26 ,297.313 ],
[ -1098.6 ,296.579 ],
[ -1098.71 ,335.984 ],
[ -1097.5 ,343.588 ],
[ -1096.84 ,364.1 ],
[ -1095.99 ,382.726 ],
[ -1095.04 ,395.717 ],
[ -1094.19 ,412.209 ],
[ -1093.07 ,422.79 ],
[ -1092.33 ,441.698 ],
[ -1091.47 ,458.251 ],
[ -1090.59 ,474.092 ],
[ -1089.83 ,493.033 ],
[ -1088.75 ,504.253 ],
[ -1088.09 ,525.595 ],
[ -1087.41 ,545.61 ],
[ -1086.48 ,561.185 ],
[ -1085.72 ,579.391 ],
[ -1084.83 ,594.811 ],
[ -1083.23 ,593.258 ],
[ -1082.5 ,613.172 ],
[ -1081.39 ,623.078 ],
[ -1079.95 ,625.454 ],
[ -1079.62 ,655.043 ],
[ -1078.17 ,656.402 ],
[ -1076.86 ,661.297 ],
[ -1075.61 ,668.208 ],
[ -1075.28 ,697.284 ],
[ -1073.34 ,686.761 ],
[ -1072.24 ,696.537 ],
[ -1070.2 ,684.011 ],
[ -1069.16 ,694.997 ],
[ -1067.38 ,688.451 ],
[ -1066.05 ,693.354 ],
[ -1064.47 ,691.191 ],
[ -1063.78 ,712.209 ],
[ -1062.57 ,719.87 ],
[ -1061.66 ,734.032 ],
[ -1060.42 ,741.598 ],
[ -1059.67 ,760.531 ],
[ -1058.33 ,764.734 ],
[ -1057.18 ,774.286 ],
[ -1055.49 ,769.233 ],
[ -1055.97 ,818.704 ],
[ -1054.82 ,827.687 ],
[ -1053.02 ,820.618 ],
[ -1052.44 ,843.83 ],
[ -1051.12 ,848.895 ],
[ -1050.86 ,879.231 ],
[ -1050.84 ,916.502 ],
[ -1049.51 ,921.361 ],
[ -1047.27 ,902.534 ],
[ -1046.63 ,924.385 ],
[ -1045.74 ,939.488 ],
[ -1044.5 ,946.798 ],
[ -1043.81 ,967.127 ],
[ -1042.08 ,961.49 ],
[ -1041.74 ,990.59 ],
[ -1041.12 ,1013.1 ],
[ -1039.39 ,1007.66 ],







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
plt.savefig(save_path, bbox_inches='tight')
plt.show()
