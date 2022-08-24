from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate
from scipy.interpolate import make_interp_spline, BSpline
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
x_axis = np.array(Extract(list_val,0))
y_axis = np.array(Extract(list_val,1))
xnew = np.linspace(x_axis.min(), x_axis.max(), 300) 

spl = make_interp_spline(x_axis, y_axis, k=2)  # type: BSpline
power_smooth = spl(xnew)

plt.plot(xnew, power_smooth, color='brown')

plt.scatter(x_axis,y_axis,color='black')
plt.xlabel('Cutoff radius in Å')
plt.ylabel('Time in fs')
plt.title('Cutoff radius vs Time in 923 Gold cluster ')
plt.grid()
# plt.show()
plt.savefig(data_file, pad_inches=1)