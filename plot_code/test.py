import numpy as np
import matplotlib.pyplot as plt
import os
# from scipy.interpolate import spline

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 3871
timestep_1 = 1 
tau = 100 * timestep_1
delQ = 250

list_val = [
[-14139.4 , 114.014],
[-13988.9 , 489.945],
[-13840.4 , 604.967],
[-13692 , 820.849],
[-13542.8 , 986.617],
[-13394.8 , 917.345],
[-13247.7 , 841.949],
[-13101.3 , 934.129],
[-12957.9 , 978.34],
[-12814.3 , 1044.38],
[-12670.7 , 1082.32],
[-12527 , 1250.52],
[-12387.7 , 1276.05],
[-12248.6 , 1643.3],
[-12106.1 , 1756.4],
[-11967.8 , 1733.62],
[-11832.1 , 1968.12],
[-11693.9 , 2222.4],
[-11563.6 , 2391.34],
[-11435.8 , 2323.2],

]

#region plot start
x_axis = Extract(list_val,1)
y_axis = Extract(list_val,0)



plt.xlabel("Temperature (K)")
plt.ylabel("Energy (eV)")
plt.suptitle(f"Energy vs Temperature ( {atoms_num} atoms )")
plt.title(f"Tau = {tau} fs, timestep = {timestep_1} fs, ΔQ = {delQ} eV")
plt.legend(f"{timestep_1} fs", loc = 2)
plt.plot(x_axis,y_axis)
plt.grid()
#endregion

# Save fig

path = os.path.join(f"{atoms_num} atoms")
save_path = os.path.join(path,f"{tau}_{timestep_1}_{delQ}_EvsT.png")
os.makedirs(path, exist_ok=True)
# plt.savefig(save_path, bbox_inches='tight')
plt.show()
