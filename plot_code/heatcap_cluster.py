from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import make_interp_spline, BSpline
def Extract(lst, index):
    return [item[index] for item in lst]

project_path = os. getcwd()

list_val = [
[ 13 ,0.00190055 ],
[ 55 ,0.0128615 ],
[ 147 ,0.0384226 ],
[ 309 ,0.111419 ],
[ 923 ,0.196057 ],
[ 2869 ,1.71892 ],
[ 3871 ,0.882341 ],
[ 5083 ,1.49819 ],
[ 8217 ,1.37581 ],
[ 10179 ,1.48248 ],
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
plt.ylabel(r'Heat capacity $C_p$ (eV/K)')
plt.title('Heat capacity vs cluster size')
plt.grid(True,which="both")
path = "plot_code/milestone_plots/"
save_path = os.path.join(path,f"heatcap_cluster_lin.png")
plt.legend(['Curvefit','Real data'])
os.makedirs(path, exist_ok=True)
# plt.show()
plt.savefig(save_path, pad_inches=1)
