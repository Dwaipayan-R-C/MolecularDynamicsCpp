from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from scipy.interpolate import make_interp_spline
from scipy.ndimage.filters import gaussian_filter1d

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 3871
timestep_1 = 2 
tau = 600 * timestep_1
delQ = 230 #eV


list_val = [
#region start
[ -14170.8 ,118.201 ],
[ -14125.9 ,399.557 ],
[ -13951.8 ,389.923 ],
[ -13895.2 ,510.078 ],
[ -13889.6 ,564.1 ],
[ -13950.5 ,748.04 ],
[ -13802.6 ,568.406 ],
[ -13878.6 ,710.564 ],
[ -13817.2 ,709.528 ],
[ -13784.7 ,663.163 ],
[ -13829.6 ,757.808 ],
[ -13783.4 ,712.008 ],
[ -13778 ,714.7 ],
[ -13766.9 ,738.988 ],
[ -13741.2 ,718.026 ],
[ -13744.5 ,744.748 ],
[ -13718.7 ,736.862 ],
[ -13713.6 ,741.798 ],
[ -13703.1 ,750.582 ],
[ -13684.7 ,736.921 ],
[ -13681.2 ,753.847 ],
[ -13664.1 ,747.041 ],
[ -13660.5 ,752.237 ],
[ -13644.3 ,744.929 ],
[ -13640.8 ,756.128 ],
[ -13624.6 ,745.591 ],
[ -13619.4 ,748.213 ],
[ -13606.6 ,744.335 ],
[ -13604.8 ,754.141 ],
[ -13592.2 ,753.78 ],
[ -13589.6 ,760.29 ],
[ -13579 ,755.265 ],
[ -13572.6 ,755.83 ],
[ -13568.2 ,761.354 ],
[ -13558.4 ,753.866 ],
[ -13548.9 ,749.414 ],
[ -13549.1 ,761.487 ],
[ -13541.5 ,760.753 ],
[ -13535.6 ,760.224 ],
[ -13530.7 ,762.474 ],
[ -13527.2 ,768.479 ],
[ -13517.1 ,761.57 ],
[ -13519.8 ,771.331 ],
[ -13515.7 ,771.043 ],
[ -13519.1 ,781.854 ],
[ -13513 ,775.687 ],
[ -13513.1 ,781.676 ],
[ -13511.6 ,781.971 ],
[ -13510.4 ,785.115 ],
[ -13512.3 ,795.213 ],
[ -13511 ,802.908 ],
[ -13517.4 ,825.397 ],
[ -13520.7 ,843.974 ],
[ -13533.2 ,886.367 ],
[ -13544 ,944.805 ],
[ -13557.4 ,1017.28 ],
[ -13574.5 ,1117.46 ],
[ -13596.1 ,1248.74 ],
[ -13624.5 ,1414.89 ],
[ -13638 ,1582.32 ],
[ -13649.8 ,1760.15 ],
[ -13641.8 ,1916.62 ],
[ -13619.2 ,2056.06 ],
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
plt.ylabel("Potential Energy (eV)")
plt.suptitle(f"Potential Energy vs Temperature ( {atoms_num} atoms )")
plt.title(f"Tau = {tau} fs, timestep = {timestep_1} fs, ΔQ = {delQ} eV")
# plt.plot(x_axis,y_axis, color='crimson')
# plt.plot(x_axis,f2)

plt.grid()
#endregion

# Save fig

path = os.path.join(f"plot_code/{atoms_num} atoms")
save_path = os.path.join(path,f"{tau}_{timestep_1}_{delQ}_PotentialEnergyVsTemperature.png")
os.makedirs(path, exist_ok=True)
# plt.text(x_axis[0],-3100,f'Heat Capacity - {heat_cap} eV/K \nMelting point - 750 K \nLatent heat - 127.5 eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )
# plt.xticks(np.arange(min(x_axis), max(x_axis)+1, 200))
plt.savefig(save_path, bbox_inches='tight')
plt.show()
