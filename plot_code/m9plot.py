from turtle import color
import matplotlib.pyplot as plt
import csv
import os

i = []
s = []
r = []
l = []
cuurent_dir = os.getcwd()
list_dir = ["data\\milestone9.dat",'data\\milestone9_50.dat','data\\milestone9_0_large.dat']
save_path = os.path.join(cuurent_dir,list_dir[2])

with open(save_path, 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    next(plots)
    for row in plots:
        i.append(float(row[0]))
        s.append(float(row[1]))
        r.append(float(row[2]))
        l.append(float(row[3]))
        
total_force=[]

for v in range(0,len(r)):
    total_force.append((-r[v]+l[v])/2)

strain_value = []
force_strain=[]
for j in range(100,len(s),200):
    strain_value.append(s[j])
    force_strain.append(total_force[j])
   
fig, axes = plt.subplots(1, 2, figsize=(11,4))
ax1 = axes[0]
ax2=axes[1]

default_x_ticks = range(len(i))
ax1.plot(i,total_force)
ax2.plot(strain_value,force_strain, 'r-o')
















ax1.legend()
ax2.legend()
ax1.grid()
ax2.set_xlabel('Strain')
ax1.set_xlabel('Timestep (fs)')
ax1.set_ylabel('Force (eV/Å)')
ax2.set_ylabel('Force (eV/Å)')
ax2.grid()
ax2.set_title('Force vs strain behavior (Small whisker)')
ax1.set_title('Force vs timestep (Small whisker)')
ax1.legend(['Temp 0K'])
ax2.legend(['Temp 0K'])
path = "plot_code/milestone_plots/"
save_path = os.path.join(path,f"m9_temp0_plot_small.png")
plt.savefig(save_path)
plt.tight_layout()
plt.show()