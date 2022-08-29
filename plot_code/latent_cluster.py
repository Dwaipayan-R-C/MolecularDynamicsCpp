import numpy as np
import matplotlib.pyplot as plt
import os
def Extract(lst, index):
    return [item[index] for item in lst]

project_path = os. getcwd()
list_val = [
[ 147 ,6.07],
[ 309 , 16.17],
[ 923 ,47.69 ],
[ 1415 ,105.9 ],
[ 2057 ,105.86 ],
[ 2869 ,151.42 ],
[ 3871 ,212.1],
[ 5083 ,484.1 ],
[ 8217 ,604.7],
[ 10179 ,766.3],
]

#region plot start
x_axis = Extract(list_val,0)
y_axis = Extract(list_val,1)
default_x_ticks = range(len(x_axis))
plt.plot(default_x_ticks, y_axis, color='brown')
plt.xticks(default_x_ticks, x_axis)
plt.scatter(default_x_ticks, y_axis, color='black')
plt.xlabel('Cluster size')
plt.ylabel(r'Latent heat (eV)')
plt.title('Latent heat vs cluster size')
plt.grid(True,which="both")
path = "plot_code/milestone_plots/"
save_path = os.path.join(path,f"latent_cluster_lin.png")

os.makedirs(path, exist_ok=True)
plt.show()
plt.savefig(save_path, pad_inches=1)
