from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from scipy.interpolate import make_interp_spline
from scipy.ndimage.filters import gaussian_filter1d

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 8217
timestep_1 = 2 
tau = 1000 * timestep_1
delQ = 150 #eV


list_val = [
#region start
[ -30024.4 ,191.088 ],
[ -29869.4 ,218.227 ],
[ -29721.8 ,318.132 ],
[ -29577.8 ,406.272 ],
[ -29427.6 ,487.932 ],
[ -29273.3 ,514.939 ],
[ -29122.2 ,560.689 ],
[ -28974.7 ,669.204 ],
[ -28825.4 ,741.094 ],
[ -28670.3 ,754.919 ],
[ -28518.1 ,776.489 ],
[ -28369.3 ,852.518 ],
[ -28219.3 ,913.076 ],
[ -28068 ,939.386 ],
[ -27916.1 ,935.225 ],
[ -27765 ,959.615 ],
[ -27614.6 ,990.44 ],
[ -27463.3 ,1002.13 ],
[ -27312 ,1029.7 ],
[ -27161.2 ,1077.79 ],
[ -27010 ,1125.89 ],
[ -26859 ,1174.39 ],
[ -26708.3 ,1232.12 ],
[ -26556.8 ,1276.82 ],
[ -26406.1 ,1335.33 ],




#endRegion
]

curve_list_1 = [
    [ -29122.2 ,560.689 ],
[ -28974.7 ,669.204 ],
[ -28825.4 ,741.094 ],
[ -28670.3 ,754.919 ],
[ -28518.1 ,776.489 ],
[ -28369.3 ,852.518 ],
[ -28219.3 ,913.076 ],
[ -28068 ,939.386 ],
]


curve_list_2 = [
    [ -27463.3 ,1002.13 ],
[ -27312 ,1029.7 ],
[ -27161.2 ,1077.79 ],
[ -27010 ,1125.89 ],
[ -26859 ,1174.39 ],
[ -26708.3 ,1232.12 ],
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
plt.plot(x_axis,y_axis, color = 'brown')
plt.scatter(x_axis,y_axis,  color = 'black')
# plt.plot(curve_fit_x, a*curve_fit_x+b, color='blue', linewidth=1.5)
# plt.plot(curve_fit_x_high, a_high*curve_fit_x_high+b_high, color='blue', linewidth=1.5)
plt.grid()
#endregion

plt.text(x_axis[0],-27500,f'Heat Capacity - {round(a,5)} eV/K \nMelting point - {melting} K \nLatent heat - {round((28068-27463.3),2)} eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )
plt.legend(["Simulated","Curvefit"])
path = os.path.join(f"plot_code/cluster_sizes")
os.makedirs(path, exist_ok=True)
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
# plt.savefig(save_path, bbox_inches='tight')
plt.show()


