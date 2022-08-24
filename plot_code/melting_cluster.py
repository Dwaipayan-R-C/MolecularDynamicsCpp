import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import make_interp_spline, BSpline
def Extract(lst, index):
    return [item[index] for item in lst]

project_path = os. getcwd()

list_val = [[ 13 ,320 ],
[ 55 ,470 ],
[ 147 ,630 ],
[ 309 ,650 ],
[ 923 ,800 ],
[ 3871 ,830 ],
[ 5083 ,800 ],
[ 8217 ,850 ],
[ 10179 ,850 ],
]

#region plot start
x_axis = np.array(Extract(list_val,0))
y_axis = np.array(Extract(list_val,1))
xnew = np.linspace(x_axis.min(), x_axis.max(), 300) 

spl = make_interp_spline(x_axis, y_axis, k=2)  # type: BSpline
power_smooth = spl(xnew)

plt.plot(xnew, power_smooth)
plt.xscale('log')
plt.yscale('log')
plt.plot(x_axis,y_axis, 'r-o')
plt.xlabel('Cluster size')
plt.ylabel('Melting temperature (K)')
plt.title('Melting point vs cluster size')
plt.grid(True,which="both")
path = "plot_code/milestone_plots/"
save_path = os.path.join(path,f"melting_cluster.png")
plt.legend(['Curvefit','Real data'])
os.makedirs(path, exist_ok=True)
# plt.show()
plt.savefig(save_path, pad_inches=1)
