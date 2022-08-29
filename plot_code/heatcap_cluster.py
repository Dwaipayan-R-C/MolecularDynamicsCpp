
import numpy as np
import matplotlib.pyplot as plt
import os
def Extract(lst, index):
    return [item[index] for item in lst]

project_path = os. getcwd()
list_val = [
[ 147 ,0.05171],
[ 309 , 0.08977],
[ 923 ,0.1762 ],
[ 1415 ,0.4633 ],
[ 2057 ,0.6721 ],
[ 2869 ,0.98359 ],
[ 3871 ,1.19949 ],
[ 5083 ,1.57704 ],
[ 8217 ,2.87633 ],
[ 10179 ,3.63513 ],
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
