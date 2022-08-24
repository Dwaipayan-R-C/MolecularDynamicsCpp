from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import make_interp_spline, BSpline
def Extract(lst, index):
    return [item[index] for item in lst]

project_path = os. getcwd()

list_val = [
[ 13 ,0.37 ],
[ 55 ,2 ],
[ 147 ,5 ],
[ 309 ,20 ],
[ 923 ,50 ],
[ 3871 ,270 ],
[ 5083 ,600 ],
[ 8217 ,900 ],
[ 10179 ,1100 ],
]

#region plot start
x_axis = np.array(Extract(list_val,0))
y_axis = np.array(Extract(list_val,1))
xnew = np.linspace(x_axis.min(), x_axis.max(), 300) 

spl = make_interp_spline(x_axis, y_axis, k=2)  # type: BSpline
power_smooth = spl(xnew)

plt.plot(xnew, power_smooth, color='brown')
# plt.xscale('log')
# plt.yscale('log')
plt.scatter(x_axis,y_axis, color='black')
plt.xlabel('Cluster size')
plt.ylabel(r'Latent heat (eV/)')
plt.title('Latent heat vs cluster size')
plt.grid(True,which="both")
path = "plot_code/milestone_plots/"
save_path = os.path.join(path,f"latent_cluster_lin.png")
plt.legend(['Curvefit','Real data'])
os.makedirs(path, exist_ok=True)
# plt.show()
plt.savefig(save_path, pad_inches=1)
