from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from scipy.interpolate import make_interp_spline
from scipy.ndimage.filters import gaussian_filter1d

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 10179
timestep_1 = 2 
tau = 1000 * timestep_1
delQ = 190 #eV


list_val = [
#region start
[ -37202 ,183.08 ],
[ -37009.2 ,243.961 ],
[ -36830.1 ,354.09 ],
[ -36646.2 ,451.987 ],
[ -36446.6 ,455.181 ],
[ -36255.4 ,537.448 ],
[ -36068.7 ,622.115 ],
[ -35873.4 ,656.48 ],
[ -35680.6 ,697.632 ],
[ -35494.4 ,814.202 ],
[ -35302.4 ,854.839 ],
[ -35108.4 ,855.206 ],
[ -34918.9 ,925.177 ],
[ -34728.1 ,966.904 ],
[ -34535.7 ,962.687 ],
[ -34344.2 ,987.116 ],
[ -34153.2 ,1003.85 ],
[ -33961.8 ,1020.98 ],
[ -33770.2 ,1051.45 ],
[ -33579.5 ,1109.24 ],
[ -33387.8 ,1155.74 ],
[ -33196.8 ,1209.03 ],
[ -33005.3 ,1263.49 ],
[ -32814.2 ,1316.43 ],



#endRegion
]

curve_list_1 = [
#region start
[ -36068.7 ,622.115 ],
[ -35873.4 ,656.48 ],
[ -35680.6 ,697.632 ],
[ -35494.4 ,814.202 ],
[ -35302.4 ,854.839 ],
[ -35108.4 ,855.206 ],
[ -34918.9 ,925.177 ],
[ -34728.1 ,966.904 ],
#endRegion
]

curve_list_2 = [
#region start
[ -33961.8 ,1020.98 ],
[ -33770.2 ,1051.45 ],
[ -33579.5 ,1109.24 ],
[ -33387.8 ,1155.74 ],
[ -33196.8 ,1209.03 ],
[ -33005.3 ,1263.49 ],
[ -32814.2 ,1316.43 ],

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
# plt.legend(f"{timestep_1} fs", loc = 2)
# plt.plot(smoothed_mode(floor_range), floor_range,ls='solid', color='crimson')

plt.plot(x_axis,y_axis, color = 'brown')
plt.scatter(x_axis,y_axis,  color = 'black')
plt.plot(curve_fit_x, a*curve_fit_x+b, color='blue', linewidth=1.5)
plt.plot(curve_fit_x_high, a_high*curve_fit_x_high+b_high, color='blue', linewidth=1.5)
plt.grid()
#endregion

path = os.path.join(f"{atoms_num} atoms")
save_path = os.path.join(path,f"{tau}_{timestep_1}_{delQ}_EvsT.png")
os.makedirs(path, exist_ok=True)
plt.text(x_axis[0],-34000,f'Heat Capacity - {round(a,5)} eV/K \nMelting point - {melting} K \nLatent heat - {round((34728.1-33961.8),2)} eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )
plt.legend(["Simulated","Curvefit"])
# plt.show()
# # Save fig
path = os.path.join(f"plot_code/cluster_sizes")
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
plt.savefig(save_path, bbox_inches='tight')
# plt.show()
