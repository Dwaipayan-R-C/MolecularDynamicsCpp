import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
# from scipy.interpolate import spline

def Extract(lst, index):
    return [item[index] for item in lst]

project_path = os. getcwd()
data_file = project_path + '/data/milestone5.dat'
list_val = np.loadtxt(data_file, unpack = True)
list_val = []
with open(data_file, 'r') as a:
    for line in a:
        b = [list(line.strip())]
        list_val.append(line)
        # print(b)
    # 
list_val = [s.replace("\n", "") for s in list_val]
# list_val = list(map(list, list_val))
print(list_val)

#region plot start
x_axis = Extract(list_val,1)
y_axis = Extract(list_val,0)
# smoothed_mode = interpolate.interp1d(y_axis, x_axis, 'cubic')
# floor_range = np.linspace(min(y_axis), max(y_axis), 500)

plt.plot(x_axis,y_axis,  'r-o')
plt.grid()
#endregion
# Save fig
path = os.path.join(f"miletone_plots/milestone5")
save_path = os.path.join(path,f"berendsen.png")
os.makedirs(path, exist_ok=True)
plt.legend(["Simulated","Curvefit"])
plt.savefig(save_path, bbox_inches='tight')
plt.show()