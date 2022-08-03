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
delQ = 10 #eV


list_val = [
#region start
[ -3352.6 ,128.98 ],
[ -3346.93 ,167.277 ],
[ -3343.01 ,218.274 ],
[ -3337.95 ,262.104 ],
[ -3331.24 ,290.664 ],
[ -3325.66 ,329.388 ],
[ -3321.78 ,378.17 ],
[ -3315.1 ,406.572 ],
[ -3309.78 ,447.912 ],
[ -3304.19 ,486.743 ],
[ -3298.81 ,525.791 ],
[ -3293.14 ,562.198 ],
[ -3286.83 ,594.343 ],
[ -3280.55 ,626.2 ],
[ -3275.77 ,670.369 ],
[ -3267.91 ,689.276 ],
[ -3261 ,715.844 ],
[ -3255.21 ,751.628 ],
[ -3248.72 ,781.764 ],
[ -3241.17 ,803.058 ],
[ -3232.99 ,819.439 ],
[ -3227.37 ,856.037 ],
[ -3217.4 ,857.123 ],
[ -3209.12 ,872.246 ],
[ -3199.94 ,879.731 ],
[ -3190.23 ,883.032 ],
[ -3179.68 ,879.209 ],
[ -3173.78 ,913.602 ],
[ -3167.03 ,941.791 ],
[ -3159.71 ,964.769 ],
[ -3152.21 ,985.928 ],
[ -3144.99 ,1010.24 ],
[ -3140.27 ,1055 ],
[ -3135.02 ,1095.25 ],
[ -3128.48 ,1124.62 ],
[ -3120.28 ,1140.08 ],
[ -3113.96 ,1171.96 ],
[ -3107.25 ,1199.96 ],
[ -3101.73 ,1237.65 ],
[ -3094.52 ,1261.29 ],
[ -3090.85 ,1314.06 ],
[ -3084.9 ,1349.26 ],
[ -3076.16 ,1360.53 ],

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
# plt.savefig(save_path, bbox_inches='tight')
plt.show()