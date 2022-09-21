from turtle import color
import matplotlib.pyplot as plt
import csv
import os


cuurent_dir = os.getcwd()
list_dir = ["data\\milestone9_0_002_small.dat","data\\milestone9_0_004_small.dat","data\\milestone9_0_01_small.dat"]
# list_dir = ["data\\milestone9_0.dat",'data\\milestone9_100.dat']
temp_list=['2e8','4e8','10e8']
fig, axes = plt.subplots(1, 1)

for data in range(0,len(temp_list)):  
    i = []
    s = []
    r = []
    l = []  
    # ax1 = axes[0]
    ax2=axes
    save_path = os.path.join(cuurent_dir,list_dir[data])
    with open(save_path, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        next(plots)
        for row in plots:
            i.append(float(row[0]))
            s.append(float(row[1]))
            r.append(float(row[2]))
            l.append(float(row[3]))
        
    total_force=[]
    total=[]

    for v in range(0,len(r)):
        total_force.append((-r[v]+l[v])/2)
    for v in range(100,len(r)):
        total.append((-r[v]+l[v])/2)

    strain_value = []
    force_strain=[]
    calculated_val=0
    avg_strain=0
    avg_force=0
    for j in range(100,len(s)):
        if(calculated_val>=5 or calculated_val<10):
            avg_strain+=s[j]
            avg_force+=(total_force[j]/2500*160.21)
            if(calculated_val==9):
                strain_value.append(avg_strain/5)
                force_strain.append(avg_force/5) 
                avg_strain=0
                avg_force=0
                calculated_val=0               
        calculated_val+=1
   
    # ax1.grid()
    # ax2.grid()
    # ax1.plot(i[100:],total)
    ax2.plot(strain_value,force_strain, '-*')
    # ax1.legend([f'Temp {temp_list[data]}K'])

    ax2.set_xlabel('Strain')
    # ax1.set_xlabel('Timestep (fs)')
    # ax1.set_ylabel('Force (eV/Å)')
    ax2.set_ylabel('Stress (GPa)')
    ax2.set_title('Force vs strain behavior (Small whisker)')
    # ax1.set_title('Force vs timestep (Small whisker)')

ax2.legend([f'Strain rate {temp_list[0]}'+r' s$^{-1}$ ', f'Strain rate {temp_list[1]}'+r' s$^{-1}$ ', f'Strain rate {temp_list[2]}'+r' s$^{-1}$ '])
path = "plot_code/milestone_plots/"
save_path = os.path.join(path,f"m9_0_plot_small_strainrate.png")
ax2.grid()
plt.tight_layout()
# plt.savefig(save_path)

plt.show()