import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
# from scipy.interpolate import spline

def Extract(lst, index):
    return [item[index] for item in lst]



list_val = [
    [-3315.41 , 0.0358755 , -3315.37],
[-3343.02 , 30.2101 , -3312.81],
[-3342.98 , 31.4397 , -3311.54],
[-3340.56 , 29.7143 , -3310.84],
[-3336.79 , 26.026 , -3310.76],
[-3338.05 , 27.0876 , -3310.96],
[-3340.2 , 28.9633 , -3311.24],
[-3340.03 , 28.8556 , -3311.17],
[-3340.31 , 28.9243 , -3311.39],
[-3341.11 , 29.8139 , -3311.29],
[-3340.89 , 29.3227 , -3311.57],
[-3343.93 , 32.8448 , -3311.08],
[-3342.24 , 31.0965 , -3311.14],
[-3340.25 , 29.3049 , -3310.94],
[-3340.74 , 29.7193 , -3311.02],
[-3341.84 , 30.7983 , -3311.04],
[-3342.39 , 31.2076 , -3311.18],
[-3342.65 , 31.2875 , -3311.36],
[-3341.94 , 30.6266 , -3311.31],
[-3342.13 , 30.919 , -3311.22],
[-3342.18 , 31.0581 , -3311.13],
[-3341.69 , 30.5725 , -3311.12],
[-3342.08 , 30.9265 , -3311.15],
[-3342.09 , 30.9064 , -3311.19],
[-3342.08 , 30.8484 , -3311.23],
[-3342.85 , 31.6767 , -3311.17],
[-3341.71 , 30.553 , -3311.16],
[-3342.86 , 31.7028 , -3311.16],
[-3343.41 , 32.2627 , -3311.15],
[-3342.42 , 31.2077 , -3311.22],
[-3342.45 , 31.2926 , -3311.15],
[-3342.08 , 30.8893 , -3311.19],
[-3343.14 , 31.9357 , -3311.21],
[-3342.33 , 31.1635 , -3311.17],
[-3342.15 , 31.0916 , -3311.05],
[-3342.79 , 31.6696 , -3311.12],
[-3341.8 , 30.7348 , -3311.07],
[-3342.54 , 31.4274 , -3311.11],
[-3342.64 , 31.5607 , -3311.08],
[-3342.98 , 31.8332 , -3311.14],
[-3342.43 , 31.2954 , -3311.14],
[-3342.32 , 31.1034 , -3311.21],
[-3342.83 , 31.6314 , -3311.2],
[-3342.06 , 30.906 , -3311.15],
[-3342.79 , 31.6845 , -3311.1],
[-3342.35 , 31.2009 , -3311.15],
[-3341.32 , 30.2833 , -3311.04],
[-3341.91 , 30.8406 , -3311.07],
[-3343.03 , 31.9004 , -3311.13],
[-3343.12 , 31.8643 , -3311.25],



]


#region plot start

kinetic = Extract(list_val,0)
potenrial = Extract(list_val,1)
total = Extract(list_val,2)
x_axis = np.linspace(0,2,len(kinetic))
# smoothed_mode = interpolate.interp1d(y_axis, x_axis, 'cubic')
# floor_range = np.linspace(min(y_axis), max(y_axis), 500)

plt.plot(x_axis,kinetic,  color='green')
plt.plot(x_axis,potenrial,  color='blue')
plt.plot(x_axis,total,  color='red')
plt.grid()
#endregion
# Save fig
path = "plot_code/milestone_plots/"
save_path = os.path.join(path,f"mpi_milestone8.png")
os.makedirs(path, exist_ok=True)
plt.legend(["Kinetic Energy","Potential Energy", "Total Energy"])
plt.xlabel('Timestep (LJ Unit)')
plt.ylabel('Energy (LJ Unit)')
plt.title("NSE Energy conservation for 4 processes")
plt.savefig(save_path, bbox_inches='tight')
plt.show()