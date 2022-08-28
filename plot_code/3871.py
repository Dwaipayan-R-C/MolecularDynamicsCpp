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
[ -14164.4 ,120.832 ],
[ -14092 ,181.349 ],
[ -14024.4 ,240.76 ],
[ -13956.1 ,316.989 ],
[ -13886.1 ,393.936 ],
[ -13816.2 ,474.837 ],
[ -13746 ,548.582 ],
[ -13675.9 ,614.534 ],
[ -13604.7 ,682.283 ],
[ -13533.6 ,721.392 ],
[ -13462.3 ,720.279 ],
[ -13391.8 ,775.748 ],
[ -13321.2 ,838.524 ],
[ -13250.6 ,870.345 ],
[ -13179.9 ,900.036 ],
[ -13109.4 ,914.958 ],
[ -13038.7 ,926.147 ],
[ -12968.3 ,956.188 ],
[ -12897.8 ,977.661 ],
[ -12827.4 ,1015.65 ],
[ -12757 ,1063.15 ],
[ -12686.6 ,1108.98 ],
[ -12616.1 ,1162.23 ],
[ -12545.7 ,1222.83 ],
[ -12475.2 ,1274.05 ],
[ -12404.5 ,1311.96 ],



#endRegion
]
curve_list_1 = [
    
[ -13689.9 ,616.14 ],
[ -13618.8 ,662.734 ],
[ -13547.9 ,714.162 ],
[ -13476.7 ,758.058 ],
[ -13405.9 ,798.895 ],
[ -13335.1 ,834.991 ],
[ -13264.7 ,880.539 ],
[ -13194.1 ,909.781 ],

]


curve_list_2=[
    [-12982 ,929.085 ],
[ -12911.5 ,959.527 ],
[ -12841.2 ,1002.17 ],
[ -12770.8 ,1049.71 ],
[ -12700.4 ,1105.52 ],
[ -12630 ,1158.58 ],
[ -12559.5 ,1214 ],
[ -12489.1 ,1259.85 ],
[ -12418.5 ,1307.63 ],
[ -12348.1 ,1370.53 ],
[ -12277.7 ,1416.56 ],
[ -12207.3 ,1474.43 ],
[ -12136.4 ,1513.68 ],
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
# plt.plot(curve_fit_x, a*curve_fit_x+b, color='blue', linewidth=1.5)
# plt.plot(curve_fit_x_high, a_high*curve_fit_x_high+b_high, color='blue', linewidth=1.5)
plt.grid()
#endregion

path = os.path.join(f"{atoms_num} atoms")
save_path = os.path.join(path,f"{tau}_{timestep_1}_{delQ}_EvsT.png")
os.makedirs(path, exist_ok=True)
# plt.text(x_axis[0],-12700,f'Heat Capacity - {round(a,5)} eV/K \nMelting point - {melting} K \nLatent heat - {round((13194.1-12982),2)} eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )
plt.legend(["Simulated","Curvefit"])
plt.show()
# # Save fig
path = os.path.join(f"plot_code/cluster_sizes")
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
# plt.savefig(save_path, bbox_inches='tight')
