import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate

def Extract(lst, index):
    return [item[index] for item in lst]

project_path = os. getcwd()
data_file = project_path + '\\plot_code\\milestone_plots\\milestone6.png'

list_val = [[ 1 ,2912 ],
[ 2 ,3047 ],
[ 3 ,3436 ],
[ 4 ,4006 ],
[ 5 ,3814 ],
[ 6 ,3672 ],
[ 7 ,4468 ],
[ 8 ,4931 ],
[ 9 ,4820 ],
[ 10 ,5391 ],]

#region plot start
x_axis = Extract(list_val,0)
y_axis = Extract(list_val,1)


plt.plot(x_axis,y_axis,  'r-o')
plt.xlabel('Cutoff radius in Å')
plt.ylabel('Time in seconds')
plt.title('Cutoff radius vs Time in 923 Gold cluster (1000 timescale)')
plt.grid()
plt.savefig(data_file, pad_inches=1)