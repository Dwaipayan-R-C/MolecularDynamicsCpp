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
delQ = 120 #eV


list_val = [
#region start
[ -30222.3 ,194.012 ],
[ -30116.8 ,212.089 ],
[ -30090.1 ,297.333 ],
[ -30054.9 ,372.131 ],
[ -29999.2 ,431.332 ],
[ -29915.4 ,469.006 ],
[ -29809.9 ,486.139 ],
[ -29770.3 ,560.611 ],
[ -29749.4 ,650.948 ],
[ -29674.6 ,693.921 ],
[ -29540.6 ,685.106 ],
[ -29467.8 ,730.58 ],
[ -29417.6 ,795.384 ],
[ -29351 ,845.88 ],
[ -29262.2 ,876.163 ],
[ -29159.1 ,893.299 ],
[ -29063.3 ,917.059 ],
[ -28979.3 ,951.345 ],
[ -28868.7 ,961.101 ],
[ -28766.6 ,978.892 ],
[ -28650.5 ,983.302 ],
[ -28545.5 ,998.327 ],
[ -28445.3 ,1017.83 ],
[ -28364.2 ,1054.9 ],
[ -28279.5 ,1089.19 ],
[ -28202.8 ,1130.41 ],
[ -28119.8 ,1166.12 ],
[ -28037.5 ,1202.71 ],
[ -27983.9 ,1265.62 ],
[ -27887.1 ,1288.53 ],
[ -27824.4 ,1343.09 ],
[ -27752.4 ,1388.88 ],
[ -27657.5 ,1413.61 ],
[ -27602.5 ,1475.24 ],
[ -27509 ,1501.52 ],






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
