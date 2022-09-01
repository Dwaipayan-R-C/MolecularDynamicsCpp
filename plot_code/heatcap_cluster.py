
import numpy as np
import matplotlib.pyplot as plt
import os
def Extract(lst, index):
    return [item[index] for item in lst]

project_path = os. getcwd()
list_val = [
[ 147 ,0.000315],
[ 309 , 0.000291],
[ 923 ,0.000318 ],
[ 1415 ,0.000305 ],
[ 2057 ,0.000306 ],
[ 2869 ,0.000311 ],
[ 3871 ,0.00031],
[ 5083 ,0.000312 ],
[ 8217 ,0.000305 ],
[ 10179 ,0.000301 ],
]

#region plot start
x_axis = Extract(list_val,0)
y_axis = Extract(list_val,1)
default_x_ticks = range(len(x_axis))
plt.plot(default_x_ticks, y_axis, color='brown')
plt.xticks(default_x_ticks, x_axis)
plt.scatter(default_x_ticks, y_axis, color='black')
plt.xlabel('Cluster size')
plt.ylabel(r'Heat capacity $C_p$ (eV/K)')
plt.title('Heat capacity vs cluster size')
plt.grid(True,which="both")
path = "plot_code/milestone_plots/"
save_path = os.path.join(path,f"heatcap_cluster_lin.png")
os.makedirs(path, exist_ok=True)
# plt.show()
plt.savefig(save_path, pad_inches=1)
