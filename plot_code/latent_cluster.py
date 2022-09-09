import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import interp1d


def Extract(lst, index):
    return [item[index] for item in lst]


project_path = os. getcwd()
list_val = [
[ 147 ,0.0413],
[ 309 , 0.0523],
[ 923 ,0.0764 ],
[ 1415 ,0.0748 ],
[ 2057 ,0.0875],
[ 2869 ,0.1085],
[ 3871 ,0.0912],
[ 5083 ,0.0952],
[ 8217 ,0.092],
[ 10179 ,0.0941],
]

#region plot start
x_axis = Extract(list_val,0)
y_axis = Extract(list_val,1)
default_x_ticks = range(len(x_axis))
plt.plot(default_x_ticks, y_axis, color='brown')
plt.xticks(default_x_ticks, x_axis)
plt.scatter(default_x_ticks, y_axis, color='black')
plt.xlabel('Cluster size')
plt.ylabel(r'Latent heat (eV/atom)')
plt.title('Latent heat vs cluster size')
plt.grid(True,which="both")
path = "plot_code/milestone_plots/"
# save_path = os.path.join(path,f"latent_cluster_lin.png")

os.makedirs(path, exist_ok=True)
plt.show()
# plt.savefig(save_path, pad_inches=1)
