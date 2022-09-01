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
[ -3337.21 ,128.956 ],
[ -3326.98 ,167.255 ],
[ -3316.98 ,218.242 ],
[ -3306.69 ,261.962 ],
[ -3296.54 ,289.661 ],
[ -3286.58 ,327.951 ],
[ -3276.55 ,372.165 ],
[ -3266.41 ,412.025 ],
[ -3256.18 ,446.903 ],
[ -3246.21 ,476.333 ],
[ -3236.28 ,516.318 ],
[ -3226 ,559.041 ],
[ -3215.81 ,588.4 ],
[ -3205.84 ,629.212 ],
[ -3195.77 ,665.129 ],
[ -3185.62 ,694.523 ],
[ -3175.57 ,722.965 ],
[ -3165.5 ,751.372 ],
[ -3155.41 ,789.207 ],
[ -3145.35 ,801.836 ],
[ -3135.25 ,819.372 ],
[ -3125.19 ,836.611 ],
[ -3115.13 ,862.186 ],
[ -3105.03 ,843.057 ],
[ -3095 ,857.611 ],
[ -3084.88 ,866.202 ],
[ -3074.82 ,893.984 ],
[ -3064.73 ,888.551 ],
[ -3054.67 ,933.131 ],
[ -3044.6 ,965.436 ],
[ -3034.52 ,971.523 ],
[ -3024.53 ,1027.53 ],
[ -3014.42 ,1052.93 ],
[ -3004.31 ,1066.9 ],
[ -2994.26 ,1106.6 ],
[ -2984.14 ,1147.43 ],
[ -2974.1 ,1152.02 ],
[ -2964.02 ,1198.81 ],
[ -2954.02 ,1250.06 ],

#endRegion
]

curve_list_1 = [
[ -3337.21 ,128.956 ],
[ -3326.98 ,167.255 ],
[ -3316.98 ,218.242 ],
[ -3306.69 ,261.962 ],
[ -3296.54 ,289.661 ],
[ -3286.58 ,327.951 ],
[ -3276.55 ,372.165 ],
[ -3266.41 ,412.025 ],
[ -3256.18 ,446.903 ],
[ -3246.21 ,476.333 ],
[ -3236.28 ,516.318 ],
[ -3226 ,559.041 ],
[ -3215.81 ,588.4 ],
[ -3205.84 ,629.212 ],
[ -3195.77 ,665.129 ],
[ -3185.62 ,694.523 ],
[ -3175.57 ,722.965 ],
[ -3165.5 ,751.372 ],
[ -3155.41 ,789.207 ],
[ -3145.35 ,801.836 ],
[ -3135.25 ,819.372 ],
[ -3125.19 ,836.611 ],
[ -3115.13 ,862.186 ],

]

curve_list_2 = [
    [ -3044.6 ,965.436 ],
[ -3034.52 ,971.523 ],
[ -3024.53 ,1027.53 ],
[ -3014.42 ,1052.93 ],
[ -3004.31 ,1066.9 ],
[ -2994.26 ,1106.6 ],
[ -2984.14 ,1147.43 ],
[ -2974.1 ,1152.02 ],
[ -2964.02 ,1198.81 ],
[ -2954.02 ,1250.06 ],
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

ax.text(0.3,0.75,f'Heat Capacity - {round(a/atoms_num,6)} eV/(K*atoms) \nMelting point - {melting} K \nLatent heat - {round((3115.13 - 3044.6)/atoms_num,4)} eV/atom' ,ha='center', va='center',fontsize=10,transform=ax.transAxes,  bbox=dict(facecolor='red', alpha=0.5) )
ax.legend(["Simulated","Curvefit"])
path = os.path.join(f"plot_code/cluster_sizes")
os.makedirs(path, exist_ok=True)
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
fig.savefig(save_path, bbox_inches='tight')
plt.show()