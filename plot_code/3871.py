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
delQ = 50 #eV


list_val = [
#region start
[ -14211.1 ,194.683 ],
[ -14187.7 ,248.168 ],
[ -14165.3 ,306.971 ],
[ -14135.8 ,350 ],
[ -14107.2 ,394.099 ],
[ -14081.6 ,442.149 ],
[ -14049.3 ,479.34 ],
[ -14021.9 ,526.748 ],
[ -13991.4 ,567.374 ],
[ -13955.7 ,596.995 ],
[ -13930.6 ,647.546 ],
[ -13898.5 ,684.292 ],
[ -13868.8 ,725.631 ],
[ -13836.7 ,762.051 ],
[ -13807 ,803.447 ],
[ -13771.8 ,833.674 ],
[ -13733.7 ,858.136 ],
[ -13695.2 ,882.189 ],
[ -13655.8 ,904.365 ],
[ -13608.7 ,911.048 ],
[ -13558.4 ,911.4 ],
[ -13515.6 ,926.271 ],
[ -13469.7 ,935.366 ],
[ -13423.5 ,943.631 ],
[ -13382 ,961.312 ],
[ -13343.7 ,985.368 ],
[ -13320.2 ,1038.9 ],
[ -13287.5 ,1074.04 ],
[ -13250.6 ,1100.75 ],
[ -13214.1 ,1128.66 ],
[ -13189.5 ,1179.86 ],
[ -13161.1 ,1223.46 ],
[ -13126.9 ,1255.93 ],
[ -13093.3 ,1289.69 ],
[ -13065.8 ,1335.06 ],
[ -13035.1 ,1374.34 ],
[ -13000.8 ,1406.36 ],
[ -12961.5 ,1428.8 ],
[ -12943.6 ,1492.57 ],




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
