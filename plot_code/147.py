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
delQ = 1 #eV


list_val = [
#region start
[ -519.783 ,63.8318 ],
[ -518.709 ,90.3085 ],
[ -517.685 ,115.549 ],
[ -516.695 ,143.794 ],
[ -515.711 ,168.664 ],
[ -514.737 ,194.178 ],
[ -513.719 ,221.59 ],
[ -512.746 ,247.172 ],
[ -511.708 ,260.458 ],
[ -510.689 ,290.93 ],
[ -509.624 ,318.592 ],
[ -508.666 ,327.918 ],
[ -507.666 ,375.555 ],
[ -506.685 ,400.948 ],
[ -505.678 ,413.018 ],
[ -504.638 ,411.831 ],
[ -503.576 ,454.964 ],
[ -502.562 ,467.08 ],
[ -501.569 ,500.905 ],
[ -500.524 ,506.531 ],
[ -499.511 ,525.962 ],
[ -498.507 ,535.296 ],
[ -497.479 ,555.81 ],
[ -496.48 ,588.277 ],
[ -495.489 ,602.178 ],
[ -494.476 ,615.968 ],
[ -493.473 ,659.154 ],
[ -492.469 ,664.059 ],
[ -491.461 ,648.863 ],
[ -490.457 ,698.811 ],
[ -489.433 ,703.181 ],
[ -488.42 ,701.773 ],
[ -487.422 ,727.073 ],
[ -486.414 ,709.996 ],
[ -485.395 ,746.457 ],
[ -484.386 ,763.671 ],
[ -483.375 ,748.972 ],
[ -482.354 ,744.148 ],
[ -481.354 ,787.119 ],
[ -480.337 ,781.885 ],
[ -479.331 ,816.393 ],
[ -478.333 ,854.667 ],
[ -477.354 ,907.795 ],
[ -476.33 ,882.212 ],
[ -475.322 ,940.869 ],
[ -474.306 ,934.872 ],
[ -473.297 ,952.659 ],
[ -472.303 ,986.254 ],
[ -471.275 ,963.088 ],
[ -470.277 ,1019.48 ],
[ -469.267 ,998.592 ],
[ -468.255 ,1014.39 ],
[ -467.274 ,1092.98 ],
[ -466.246 ,1058.12 ],
[ -465.258 ,1149.45 ],
[ -464.256 ,1169.04 ],
[ -463.254 ,1176.53 ],
[ -462.243 ,1179.18 ],
[ -461.228 ,1134.16 ],


#endRegion
]

curve_list_1 = [
#region start
[ -519.783 ,63.8318 ],
[ -518.709 ,90.3085 ],
[ -517.685 ,115.549 ],
[ -516.695 ,143.794 ],
[ -515.711 ,168.664 ],
[ -514.737 ,194.178 ],
[ -513.719 ,221.59 ],
[ -512.746 ,247.172 ],
[ -511.708 ,260.458 ],
[ -510.689 ,290.93 ],
[ -509.624 ,318.592 ],
[ -508.666 ,327.918 ],
[ -507.666 ,375.555 ],
[ -506.685 ,400.948 ],
[ -505.678 ,413.018 ],
[ -504.638 ,411.831 ],
[ -503.576 ,454.964 ],
[ -502.562 ,467.08 ],
[ -501.569 ,500.905 ],
[ -500.524 ,506.531 ],
[ -499.511 ,525.962 ],
[ -498.507 ,535.296 ],
[ -497.479 ,555.81 ],
[ -496.48 ,588.277 ],
[ -495.489 ,602.178 ],
[ -494.476 ,615.968 ],
[ -493.473 ,659.154 ],
[ -492.469 ,664.059 ],
[ -491.461 ,648.863 ],
#endRegion
]

curve_list_2 = [
#region start
[ -485.395 ,746.457 ],
[ -484.386 ,763.671 ],
[ -483.375 ,748.972 ],
[ -482.354 ,744.148 ],
[ -481.354 ,787.119 ],
[ -480.337 ,781.885 ],
[ -479.331 ,816.393 ],
[ -478.333 ,854.667 ],
[ -477.354 ,907.795 ],
[ -476.33 ,882.212 ],
[ -475.322 ,940.869 ],
[ -474.306 ,934.872 ],
[ -473.297 ,952.659 ],
[ -472.303 ,986.254 ],
[ -471.275 ,963.088 ],
[ -470.277 ,1019.48 ],
[ -469.267 ,998.592 ],
[ -468.255 ,1014.39 ],
[ -467.274 ,1092.98 ],
[ -466.246 ,1058.12 ],
[ -465.258 ,1149.45 ],
[ -464.256 ,1169.04 ],
[ -463.254 ,1176.53 ],
[ -462.243 ,1179.18 ],
[ -461.228 ,1134.16 ],


#endRegion
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
ax.plot(curve_fit_x, a*curve_fit_x+b, color='blue', linewidth=1.5)
ax.plot(curve_fit_x_high, a_high*curve_fit_x_high+b_high, color='blue', linewidth=1.5)
ax.grid()
#endregion

ax.text(0.3,0.75,f'Heat Capacity - {round(a/atoms_num,6)} eV/(K*atom) \nMelting point - {melting} K \nLatent heat - {round((491.461-485.395)/atoms_num,4)} eV/atom' ,ha='center', va='center',fontsize=10,transform=ax.transAxes,  bbox=dict(facecolor='red', alpha=0.5) )
# ax.text(0.25,0.75,f'Heat Capacity - {round(a/atoms_num,5)} eV/K \nMelting point - {melting} K \nLatent heat - {round((491.461- 485.395),2)} eV' ,ha='center', va='center',fontsize=10,transform=ax.transAxes,  bbox=dict(facecolor='red', alpha=0.5) )
ax.legend(["Simulated","Curvefit"])
path = os.path.join(f"plot_code/cluster_sizes")
os.makedirs(path, exist_ok=True)
save_path = os.path.join(path,f"{atoms_num}_Energy_Temp.png")
fig.savefig(save_path, bbox_inches='tight')
plt.show()

# plt.text(x_axis[0],-480,f'Heat Capacity - {round(a/atoms_num,5)} eV/K \nMelting point - {melting} K \nLatent heat - {round((491.461- 485.395),2)} eV' ,fontsize=10, bbox=dict(facecolor='red', alpha=0.5) )

