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
delQ = 80 #eV


list_val = [
#region start
[ -18675.7 ,198.984 ],
[ -18632.5 ,253.205 ],
[ -18585.1 ,308.76 ],
[ -18535.8 ,358.683 ],
[ -18483.1 ,404.03 ],
[ -18422.2 ,434.494 ],
[ -18374.1 ,483.803 ],
[ -18329.2 ,537.833 ],
[ -18286.9 ,594.665 ],
[ -18258.6 ,673.263 ],
[ -18217.5 ,731.964 ],
[ -18170.3 ,781.537 ],
[ -18114.4 ,819.411 ],
[ -18050.7 ,845.774 ],
[ -17986.5 ,871.044 ],
[ -17918.7 ,890.422 ],
[ -17859.8 ,923.083 ],
[ -17785.9 ,933.24 ],
[ -17718 ,952.17 ],
[ -17653 ,975.85 ],
[ -17589.7 ,1002.11 ],
[ -17520.7 ,1019.72 ],
[ -17479.8 ,1079.67 ],
[ -17424.2 ,1117.64 ],
[ -17362.3 ,1146.34 ],
[ -17324.8 ,1211.42 ],
[ -17264.1 ,1241.99 ],
[ -17207.1 ,1277.9 ],
[ -17173.6 ,1349.02 ],
[ -17113.7 ,1380.98 ],
[ -17063.8 ,1427.41 ],
[ -17020.9 ,1484.59 ],
[ -16972.1 ,1532.81 ],







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
