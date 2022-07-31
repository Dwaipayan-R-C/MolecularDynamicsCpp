from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
# from scipy.interpolate import spline

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 923
timestep_1 = 1 
tau = 1000 * timestep_1
delQ = 25 #eV


list_val = [
[-3343.96 , 124.751],
[-3303.38 , 220.678],
[-3268.74 , 352.904],
[-3244.58 , 548.144],
[-3218.67 , 678.581],
[-3200.92 , 664.755],
[-3186.79 , 619.108],
[-3168.35 , 745.226],
[-3147.67 , 826.193],
[-3128.71 , 936.977],



]

#region plot start
x_axis = np.array(Extract(list_val,1))
y_axis = np.array(Extract(list_val,0))
# smoothed_mode = interpolate.interp1d(y_axis, x_axis, 'cubic')
# floor_range = np.linspace(min(y_axis), max(y_axis), 500)

a, b = np.polyfit(x_axis, y_axis, 1)

plt.xlabel("Temperature (K)")
plt.ylabel("Energy (eV)")
plt.suptitle(f"Energy vs Temperature ( {atoms_num} atoms )")
plt.title(f"Tau = {tau} fs, timestep = {timestep_1} fs, ΔQ = {delQ} eV")
# plt.legend(f"{timestep_1} fs", loc = 2)
plt.plot(x_axis, a*x_axis+b, color='steelblue', linestyle='--', linewidth=2)
# plt.plot(a,b)
plt.grid()
#endregion

heat_cap = round((y_axis[9]-y_axis[4])/(x_axis[9]-x_axis[4]),5)
# Save fig

path = os.path.join(f"{atoms_num} atoms")
save_path = os.path.join(path,f"{tau}_{timestep_1}_{delQ}_EvsT.png")
os.makedirs(path, exist_ok=True)
# plt.text(x_axis[0],-3100,f'Heat Capacity - {heat_cap} eV/K \nMelting point - 750 K \nLatent heat - 127.5 eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )
# plt.xticks(np.arange(min(x_axis), max(x_axis)+1, 200))
# plt.savefig(save_path, bbox_inches='tight')
plt.show()
