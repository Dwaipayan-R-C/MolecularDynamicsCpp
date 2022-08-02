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
delQ = 120 #eV


list_val = [
#region start
[ -37442.9 ,183.08 ],
[ -37375.9 ,224.772 ],
[ -37374.8 ,306.14 ],
[ -37333.2 ,361.58 ],
[ -37179.5 ,340.769 ],
[ -37183.4 ,434.714 ],
[ -37138.3 ,489.141 ],
[ -37034.2 ,507.186 ],
[ -36906.5 ,505.943 ],
[ -36900.9 ,589.586 ],
[ -36876.3 ,658.633 ],
[ -36762 ,666.783 ],
[ -36652.9 ,677.574 ],
[ -36614.4 ,738.176 ],
[ -36553.8 ,783.536 ],
[ -36461.3 ,805.75 ],
[ -36380.2 ,835.785 ],
[ -36296.3 ,864.054 ],
[ -36216.5 ,894.978 ],
[ -36130.4 ,921.141 ],
[ -36020.8 ,930.029 ],
[ -35919.4 ,944.739 ],
[ -35820.2 ,961.1 ],
[ -35697.4 ,959.828 ],
[ -35584.3 ,965.806 ],
[ -35481.3 ,979.1 ],
[ -35364.4 ,982.37 ],
[ -35279 ,1009.03 ],
[ -35203.3 ,1043.35 ],
[ -35110.2 ,1064.89 ],
[ -35038.5 ,1101.78 ],
[ -34960 ,1133.89 ],
[ -34880.5 ,1165.35 ],
[ -34809.5 ,1203.02 ],
[ -34726.8 ,1232.17 ],
[ -34658 ,1271.58 ],
[ -34565.7 ,1293.42 ],
[ -34510 ,1342.75 ],
[ -34429.6 ,1373.47 ],
[ -34351.2 ,1405.66 ],
[ -34290.4 ,1451.14 ],
[ -34206 ,1479.09 ],
[ -34125 ,1509.42 ],






#endRegion
]

#region plot start
x_axis = Extract(list_val,1)
y_axis = Extract(list_val,0)
param = np.linspace(0, 1, len(x_axis))
spl = make_interp_spline(param, np.c_[x_axis,y_axis], k=2) #(1)
xnew, y_smooth = spl(np.linspace(0, 1, len(x_axis) * 100)).T #(2)
plt.plot(xnew, y_smooth, color='red')
plt.scatter(x_axis, y_axis, color="black")
plt.xlabel("Temperature (K)")
plt.ylabel("Total Energy (eV)")
plt.suptitle(f"Total Energy vs Temperature ( {atoms_num} atoms )")
plt.title(f"Tau = {tau} fs, timestep = {timestep_1} fs, ΔQ = {delQ} eV")
# plt.plot(x_axis,y_axis, color='crimson')
# plt.plot(x_axis,f2)

plt.grid()
#endregion

# Save fig

path = os.path.join(f"plot_code/cluster_sizes")
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
os.makedirs(path, exist_ok=True)
# plt.text(x_axis[0],-3100,f'Heat Capacity - {heat_cap} eV/K \nMelting point - 750 K \nLatent heat - 127.5 eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )
# plt.xticks(np.arange(min(x_axis), max(x_axis)+1, 200))
# plt.savefig(save_path, bbox_inches='tight')
plt.show()
