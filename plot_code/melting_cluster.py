import numpy as np
import matplotlib.pyplot as plt
import os
def Extract(lst, index):
    return [item[index] for item in lst]

project_path = os. getcwd()
list_val = [
[ 147 ,648],
[ 309 , 702.247],
[ 923 ,819.439],
[ 1415 ,846.884 ],
[ 2057 ,866.268 ],
[ 2869 ,893.881 ],
[ 3871 ,880.874],
[ 5083 ,884.727],
[ 8217 ,939.386],
[ 10179 ,954.244],
]

#region plot start
x_axis = Extract(list_val,0)
y_axis = Extract(list_val,1)
default_x_ticks = range(len(x_axis))
plt.plot(default_x_ticks, y_axis, color='brown')
plt.xticks(default_x_ticks, x_axis)
plt.scatter(default_x_ticks, y_axis, color='black')
plt.xlabel('Cluster size')
plt.ylabel('Melting temperature (K)')
plt.title('Melting point vs cluster size')
plt.grid(True,which="both")
path = "plot_code/milestone_plots/"
save_path = os.path.join(path,f"melting_cluster_curve.png")
os.makedirs(path, exist_ok=True)
# plt.show()
plt.savefig(save_path, pad_inches=1)
