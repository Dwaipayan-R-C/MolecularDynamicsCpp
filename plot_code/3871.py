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
tau = 1000 * timestep_1
delQ = 70 #eV


list_val = [
#region start
[ -14113.7 ,194.664 ],
[ -14043.5 ,269.779 ],
[ -13971.8 ,347.677 ],
[ -13900.7 ,417.998 ],
[ -13830.2 ,486.322 ],
[ -13760.3 ,554.13 ],
[ -13689.9 ,621.067 ],
[ -13619.3 ,684.613 ],
[ -13547.6 ,737.21 ],
[ -13476.2 ,767.625 ],
[ -13405.6 ,811.983 ],
[ -13334.8 ,847.519 ],
[ -13264 ,880.874 ],
[ -13193.4 ,910.031 ],
[ -13122.8 ,919.785 ],
[ -13052.1 ,928.272 ],
[ -12981.4 ,941.544 ],
[ -12911 ,964.48 ],
[ -12840.7 ,1005.55 ],
[ -12770.4 ,1055.83 ],
[ -12699.8 ,1101.22 ],
[ -12629.5 ,1162.68 ],
[ -12558.8 ,1204.35 ],
[ -12488.4 ,1261.17 ],

#endRegion
]
curve_list_1 = [
[ -14113.7 ,194.664 ],
[ -14043.5 ,269.779 ],
[ -13971.8 ,347.677 ],
[ -13900.7 ,417.998 ],
[ -13830.2 ,486.322 ],
[ -13760.3 ,554.13 ],
[ -13689.9 ,621.067 ],
[ -13619.3 ,684.613 ],
[ -13547.6 ,737.21 ],
[ -13476.2 ,767.625 ],
[ -13405.6 ,811.983 ],
[ -13334.8 ,847.519 ],
[ -13264 ,880.874 ],

]


curve_list_2=[
    [ -12911 ,964.48 ],
[ -12840.7 ,1005.55 ],
[ -12770.4 ,1055.83 ],
[ -12699.8 ,1101.22 ],
[ -12629.5 ,1162.68 ],
[ -12558.8 ,1204.35 ],
[ -12488.4 ,1261.17 ],
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

fig,ax = plt.subplots()
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Total Energy (eV)")
fig.suptitle(f"Total Energy vs Temperature ( {atoms_num} atoms )")
ax.set_title(f"Tau = {tau} fs, timestep = {timestep_1} fs, ΔQ = {delQ} eV")
ax.plot(x_axis,y_axis, color = 'brown')
ax.scatter(x_axis,y_axis,  color = 'black')
# ax.plot(curve_fit_x, a*curve_fit_x+b, color='blue', linewidth=1.5)
# ax.plot(curve_fit_x_high, a_high*curve_fit_x_high+b_high, color='blue', linewidth=1.5)
ax.grid()
#endregion

ax.text(0.3,0.75,f'Heat Capacity - {round(a/atoms_num,6)} eV/(K*atoms) \nMelting point - {melting} K \nLatent heat - {round((13264-12911)/atoms_num,4)} eV/atom' ,ha='center', va='center',fontsize=10,transform=ax.transAxes,  bbox=dict(facecolor='red', alpha=0.5) )
ax.legend(["Simulated","Curvefit"])
path = os.path.join(f"plot_code/cluster_sizes")
os.makedirs(path, exist_ok=True)
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
fig.savefig(save_path, bbox_inches='tight')
plt.show()
# plt.text(x_axis[0],-13000,f'Heat Capacity - {round(a/atoms_num,5)} eV/K \nMelting point - {melting} K \nLatent heat - {round((13194.1-12981.4),2)} eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )
