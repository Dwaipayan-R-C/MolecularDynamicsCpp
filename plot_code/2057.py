from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from scipy.interpolate import make_interp_spline
from scipy.ndimage.filters import gaussian_filter1d

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 2057
timestep_1 = 2 
tau = 1000 * timestep_1
delQ = 35 #eV

list_val = [
#region start
[ -7477.63 ,174.637 ],
[ -7439.07 ,205.403 ],
[ -7403.92 ,274.017 ],
[ -7370.19 ,369.222 ],
[ -7334.01 ,407.815 ],
[ -7298.63 ,463.911 ],
[ -7264.02 ,530.668 ],
[ -7228.49 ,596.663 ],
[ -7193.18 ,647.637 ],
[ -7157.69 ,687.659 ],
[ -7122.7 ,753.74 ],
[ -7087.56 ,799.627 ],
[ -7052.2 ,833.577 ],
[ -7016.87 ,866.268 ],
[ -6981.68 ,904.65 ],
[ -6946.37 ,902.478 ],
[ -6911.07 ,919.985 ],
[ -6875.82 ,933.4 ],
[ -6840.54 ,961.042 ],
[ -6805.29 ,1004.01 ],
[ -6769.98 ,1040.46 ],
[ -6734.74 ,1093.63 ],
[ -6699.46 ,1131.99 ],
[ -6664.22 ,1174.59 ],

#endRegion
]


curve_list_1 = [
#region start
[ -7477.63 ,174.637 ],
[ -7439.07 ,205.403 ],
[ -7403.92 ,274.017 ],
[ -7370.19 ,369.222 ],
[ -7334.01 ,407.815 ],
[ -7298.63 ,463.911 ],
[ -7264.02 ,530.668 ],
[ -7228.49 ,596.663 ],
[ -7193.18 ,647.637 ],
[ -7157.69 ,687.659 ],
[ -7122.7 ,753.74 ],
[ -7087.56 ,799.627 ],
[ -7052.2 ,833.577 ],
[ -7016.87 ,866.268 ],



#endRegion
]


curve_list_2 = [
#region start
[ -6840.54 ,961.042 ],
[ -6805.29 ,1004.01 ],
[ -6769.98 ,1040.46 ],
[ -6734.74 ,1093.63 ],
[ -6699.46 ,1131.99 ],
[ -6664.22 ,1174.59 ],
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

ax.text(0.3,0.75,f'Heat Capacity - {round(a/atoms_num,6)} eV/(K*atoms) \nMelting point - {melting} K \nLatent heat - {round((7016.87 - 6840.54)/atoms_num,4)} eV/atom' ,ha='center', va='center',fontsize=10,transform=ax.transAxes,  bbox=dict(facecolor='red', alpha=0.5) )
ax.legend(["Simulated","Curvefit"])
path = os.path.join(f"plot_code/cluster_sizes")
os.makedirs(path, exist_ok=True)
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
fig.savefig(save_path, bbox_inches='tight')
plt.show()

# plt.text(x_axis[0],-6900,f'Heat Capacity - {round(a/atoms_num,5)} eV/K \nMelting point - {melting} K \nLatent heat - {round((6981.68 - 6875.82),2)} eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )

