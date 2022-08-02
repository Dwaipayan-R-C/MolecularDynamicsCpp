from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from scipy.interpolate import make_interp_spline
from scipy.ndimage.filters import gaussian_filter1d

def Extract(lst, index):
    return [item[index] for item in lst]

atoms_num = 147
timestep_1 = 2 
tau = 1000 * timestep_1
delQ = 0.5 #eV


list_val = [
#region start
[ -520.981 ,64.6287 ],
[ -520.66 ,77.7725 ],
[ -520.411 ,90.6153 ],
[ -520.169 ,103.592 ],
[ -519.897 ,115.928 ],
[ -519.709 ,129.575 ],
[ -519.43 ,144.305 ],
[ -519.096 ,155.865 ],
[ -518.827 ,167.817 ],
[ -518.654 ,185.907 ],
[ -518.51 ,203.097 ],
[ -518.147 ,210.515 ],
[ -517.7 ,214.116 ],
[ -517.495 ,230.476 ],
[ -517.331 ,247.419 ],
[ -517.17 ,263.41 ],
[ -516.652 ,267.905 ],
[ -516.625 ,290.843 ],
[ -515.545 ,259.119 ],
[ -516.104 ,317.026 ],
[ -515.442 ,308.743 ],
[ -515.398 ,332.406 ],
[ -514.981 ,338 ],
[ -514.787 ,354.092 ],
[ -514.452 ,362.871 ],
[ -514.29 ,381.596 ],
[ -513.946 ,390.549 ],
[ -513.647 ,401.23 ],
[ -513.36 ,412.971 ],
[ -513.198 ,430.674 ],
[ -512.669 ,429.761 ],
[ -512.345 ,439.164 ],
[ -512.349 ,465.846 ],
[ -512.055 ,476.729 ],
[ -511.544 ,476.619 ],
[ -511.124 ,481.73 ],
[ -510.965 ,499.154 ],
[ -510.524 ,503.428 ],
[ -510.268 ,516.767 ],
[ -509.865 ,522.291 ],
[ -509.916 ,550.969 ],
[ -509.435 ,552.003 ],
[ -508.982 ,555.261 ],
[ -508.465 ,554.291 ],
[ -508.383 ,576.002 ],
[ -508.05 ,585.349 ],
[ -507.646 ,590.815 ],
[ -507.405 ,604.512 ],
[ -506.263 ,571.889 ],
[ -506.227 ,596.364 ],
[ -506.29 ,625.344 ],
[ -506.005 ,636.877 ],
[ -505.152 ,618.717 ],
[ -505.172 ,646.176 ],
[ -504.958 ,661.962 ],
[ -504.206 ,649.003 ],
[ -503.773 ,652.699 ],
[ -503.276 ,653.351 ],
[ -502.746 ,651.026 ],
[ -502.125 ,645.391 ],
[ -501.806 ,655.346 ],
[ -501.247 ,652.495 ],
[ -500.817 ,657.189 ],
[ -500.54 ,668.491 ],
[ -500.412 ,688.939 ],
[ -499.524 ,669.014 ],
[ -499.611 ,700.221 ],
[ -499.338 ,711.407 ],
[ -498.87 ,713.81 ],
[ -498.406 ,716.443 ],
[ -498.266 ,736.184 ],
[ -498.068 ,751.189 ],
[ -497.896 ,769.734 ],
[ -497.059 ,751.026 ],
[ -496.984 ,773.598 ],
[ -496.311 ,764.639 ],
[ -495.959 ,772.464 ],
[ -496.047 ,803.023 ],
[ -495.773 ,814.899 ],
[ -495.149 ,808.403 ],
[ -494.968 ,824.671 ],
[ -494.687 ,837.531 ],
[ -494.306 ,843.875 ],
[ -494.412 ,875.594 ],






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
plt.savefig(save_path, bbox_inches='tight')
plt.show()
