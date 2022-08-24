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
[ -6.70185 ,128.98 ],
[ -6.68622 ,167.277 ],
[ -6.64962 ,218.274 ],
[ -6.64899 ,262.104 ],
[ -6.664 ,290.664 ],
[ -6.64677 ,329.388 ],
[ -6.61957 ,378.17 ],
[ -6.61877 ,406.572 ],
[ -6.6146 ,447.912 ],
[ -6.58351 ,486.743 ],
[ -6.57953 ,525.782 ],
[ -6.57601 ,562.694 ],
[ -6.55674 ,596.61 ],
[ -6.54414 ,629.222 ],
[ -6.53938 ,665.819 ],
[ -6.52785 ,702.912 ],
[ -6.50834 ,724.765 ],
[ -6.50057 ,753.374 ],
[ -6.48616 ,785.862 ],
[ -6.46462 ,801.359 ],
[ -6.44197 ,801.024 ],
[ -6.42725 ,804.389 ],
[ -6.41525 ,836.796 ],
[ -6.40019 ,853.217 ],
[ -6.37853 ,866.789 ],
[ -6.35569 ,863.484 ],
[ -6.35059 ,880.757 ],
[ -6.32808 ,899.918 ],
[ -6.31667 ,914.112 ],
[ -6.30938 ,958.834 ],
[ -6.29192 ,1003.31 ],
[ -6.2882 ,1021.57 ],
[ -6.26427 ,1056.79 ],
[ -6.25157 ,1096.83 ],
[ -6.24924 ,1114.95 ],
[ -6.2384 ,1154.08 ],
[ -6.218 ,1174.42 ],
[ -6.20281 ,1209.35 ],
[ -6.18745 ,1224.49 ],
[ -6.17378 ,1246.39 ],
[ -6.15527 ,1262.57 ],
[ -6.1289 ,1341.88 ],
[ -6.14091 ,1342.68 ],
[ -6.11807 ,1401.42 ],
[ -6.11818 ,1436.93 ],



#endRegion
]

#region plot start
x_axis = Extract(list_val,1)
y_axis = Extract(list_val,0)
param = np.linspace(0, 1, len(x_axis))
spl = make_interp_spline(param, np.c_[x_axis,y_axis], k=2) #(1)
xnew, y_smooth = spl(np.linspace(0, 1, len(x_axis) * 100)).T #(2)
# plt.plot(xnew, y_smooth, color='red')
plt.plot(x_axis, y_axis,'r-o')
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
