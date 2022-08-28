from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from scipy.interpolate import make_interp_spline
from scipy.ndimage.filters import gaussian_filter1d

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 8271
timestep_1 = 2 
tau = 1000 * timestep_1
delQ = 150 #eV


list_val = [
#region start
[ -30016.3 ,194.012 ],
[ -29861.3 ,220.755 ],
[ -29713.7 ,321.619 ],
[ -29569.3 ,408.56 ],
[ -29418.4 ,477.477 ],
[ -29264.2 ,514.811 ],
[ -29113.9 ,572.099 ],
[ -28967.2 ,675.902 ],
[ -28816.3 ,744.284 ],
[ -28662.8 ,753.52 ],
[ -28510.9 ,786.257 ],
[ -28360.6 ,849.019 ],
[ -28210.7 ,911.154 ],
[ -28059.6 ,942.913 ],
[ -27907.8 ,949.517 ],
[ -27756.8 ,967.115 ],
[ -27606.5 ,1004.48 ],
[ -27454.6 ,1006.27 ],
[ -27303.8 ,1039.11 ],
[ -27152.9 ,1081.23 ],
[ -27001.5 ,1121.16 ],
[ -26851.1 ,1191.06 ],
[ -26699.5 ,1228.73 ],
[ -26549 ,1291.07 ],
[ -26397.5 ,1332.62 ],
[ -26246.8 ,1396.19 ],
[ -26095.6 ,1441.31 ],



#endRegion
]

curve_list_1 = [
    
]


curve_list_2 = [
    
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
plt.text(x_axis[0],-28000,f'Heat Capacity - {round(a,5)} eV/K \nMelting point - {melting} K \nLatent heat - {round((34728.1-33961.8),2)} eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )
plt.legend(["Simulated","Curvefit"])
plt.show()
# # Save fig
path = os.path.join(f"plot_code/cluster_sizes")
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
# plt.savefig(save_path, bbox_inches='tight')
# plt.show()

