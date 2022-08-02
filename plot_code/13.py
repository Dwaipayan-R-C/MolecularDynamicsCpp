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
delQ = 0.03 #eV


list_val = [
#region start
[ -42.198 ,208.511 ],
[ -42.1843 ,218.219 ],
[ -42.1751 ,230.598 ],
[ -42.1306 ,221.996 ],
[ -42.1526 ,252.887 ],
[ -42.1163 ,249.167 ],
[ -42.1205 ,269.542 ],
[ -42.0821 ,264.57 ],
[ -42.0669 ,273.364 ],
[ -42.0741 ,295.474 ],
[ -42.043 ,294.868 ],
[ -42.0512 ,317.586 ],
[ -42.0217 ,317.913 ],
[ -41.9908 ,317.456 ],
[ -41.9658 ,320.443 ],
[ -41.9437 ,325.185 ],
[ -41.9254 ,332.179 ],
[ -41.9003 ,335.413 ],
[ -41.8498 ,323.579 ],
[ -41.9033 ,372.677 ],
[ -41.7804 ,318.318 ],
[ -41.831 ,365.361 ],
[ -41.8474 ,392.936 ],
[ -41.8465 ,410.495 ],
[ -41.732 ,361.136 ],
[ -41.7342 ,379.926 ],
[ -41.7652 ,416.094 ],
[ -41.6131 ,347.626 ],
[ -41.5331 ,318.784 ],
[ -41.5841 ,369.355 ],
[ -41.5869 ,389.024 ],
[ -41.5274 ,369.233 ],
[ -41.4971 ,371.593 ],
[ -41.5185 ,400.558 ],
[ -41.5427 ,433.023 ],
[ -41.5404 ,449.134 ],
[ -41.495 ,441.123 ],
[ -41.447 ,429.183 ],
[ -41.4911 ,472.351 ],
[ -41.4382 ,459.135 ],
[ -41.2952 ,393.612 ],
[ -41.3751 ,458.408 ],
[ -41.4367 ,511.854 ],
[ -41.3168 ,458.059 ],
[ -41.2473 ,435.838 ],
[ -41.2574 ,458.021 ],
[ -41.2721 ,487.595 ],
[ -41.2157 ,471.538 ],
[ -41.3288 ,556.037 ],
[ -41.2821 ,545.916 ],
[ -41.1726 ,498.673 ],
[ -41.1324 ,493.72 ],
[ -41.2503 ,581.152 ],
[ -41.1282 ,530.089 ],
[ -41.2636 ,626.58 ],
[ -41.1384 ,569.16 ],
[ -41.1224 ,577.008 ],
[ -40.978 ,509.75 ],




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
