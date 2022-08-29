from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
# from scipy.interpolate import spline

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 923
timestep_1 = 2 
tau = 1000 * timestep_1
delQ = 10 #eV


list_val = [
[ -6.70185 ,128.956 ],
[ -6.68623 ,167.255 ],
[ -6.64964 ,218.242 ],
[ -6.64825 ,261.962 ],
[ -6.65634 ,289.661 ],
[ -6.65563 ,327.951 ],
[ -6.62235 ,372.165 ],
[ -6.6123 ,412.025 ],
[ -6.60994 ,446.903 ],
[ -6.60133 ,476.333 ],
[ -6.56385 ,516.318 ],
[ -6.57261 ,559.041 ],
[ -6.56277 ,588.4 ],
[ -6.53966 ,629.212 ],
[ -6.53988 ,665.129 ],
[ -6.52462 ,694.523 ],
[ -6.50595 ,722.965 ],
[ -6.49291 ,751.372 ],
[ -6.48258 ,789.207 ],
[ -6.46115 ,801.836 ],
[ -6.45694 ,819.372 ],
[ -6.44259 ,836.611 ],
[ -6.41827 ,862.186 ],
[ -6.39795 ,843.057 ],
[ -6.3729 ,857.611 ],
[ -6.3621 ,866.202 ],
[ -6.3488 ,893.984 ],
[ -6.32596 ,888.551 ],
[ -6.31951 ,933.131 ],
[ -6.30067 ,965.436 ],
[ -6.28127 ,971.523 ],
[ -6.27779 ,1027.53 ],
[ -6.26834 ,1052.93 ],
[ -6.25311 ,1066.9 ],
[ -6.23804 ,1106.6 ],
[ -6.22751 ,1147.43 ],
[ -6.21227 ,1152.02 ],
[ -6.20458 ,1198.81 ],
[ -6.19933 ,1250.06 ],
[ -6.17612 ,1275.07 ],
[ -6.15692 ,1292.83 ],
[ -6.15297 ,1324.1 ],
[ -6.12812 ,1362.27 ],
[ -6.1359 ,1395.67 ],



]


curve_list_1 = [
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


]

curve_list_2 = [

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
]

melting = np.array(Extract(curve_list_1,1))[-1]
#region plot start
x_axis = Extract(list_val,1)
y_axis = Extract(list_val,0)
# smoothed_mode = interpolate.interp1d(y_axis, x_axis, 'cubic')
# floor_range = np.linspace(min(y_axis), max(y_axis), 500)

curve_fit_x = np.array(Extract(curve_list_1,1))
curve_fit_y = np.array(Extract(curve_list_1,0))
curve_fit_x_high = np.array(Extract(curve_list_2,1))
curve_fit_y_high = np.array(Extract(curve_list_2,0))

a, b = np.polyfit(curve_fit_x, curve_fit_y, 1)
a_high, b_high = np.polyfit(curve_fit_x_high, curve_fit_y_high, 1)

plt.xlabel("Temperature (K)")
plt.ylabel("Potential Energy (eV)")
plt.suptitle(f"Potential Energy vs Temperature ( {atoms_num} atoms )")
plt.title(f"Tau = {tau} fs, timestep = {timestep_1} fs, ΔQ = {delQ} eV")
# plt.legend(f"{timestep_1} fs", loc = 2)
# plt.plot(smoothed_mode(floor_range), floor_range,ls='solid', color='crimson')

plt.plot(x_axis,y_axis,  color='brown')
plt.scatter(x_axis,y_axis,  color='black')
# plt.plot(curve_fit_x, a*curve_fit_x+b, color='black', linewidth=0.9)
# plt.plot(curve_fit_x_high, a_high*curve_fit_x_high+b_high, color='black', linewidth=0.9)
plt.grid()
#endregion


# Save fig
path = os.path.join(f"{atoms_num} atoms")
save_path = os.path.join(path,f"{tau}_{timestep_1}_{delQ}_PotentialvsTemp.png")
os.makedirs(path, exist_ok=True)
# plt.text(x_axis[0],-3150,f'Heat Capacity - {round(a,5)} eV/K \nMelting point - {melting} K \nLatent heat - {b_high - b} eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )
# plt.legend(["Simulated","Curvefit"])
# plt.xticks(np.arange(min(x_axis), max(x_axis)+1, 200))
plt.savefig(save_path, bbox_inches='tight')
# plt.show()
