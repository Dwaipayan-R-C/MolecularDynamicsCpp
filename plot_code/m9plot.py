from turtle import color
import matplotlib.pyplot as plt
import csv
import os


cuurent_dir = os.getcwd()
# list_dir = ["data\\milestone9_0.dat",'data\\milestone9_0_large.dat']
list_dir = ["data\\milestone9_0.dat",'data\\milestone9_0_large.dat','data\\milestone9_100.dat']
# temp_list=[0,0]
temp_list=[0,0,100]
fig, axes = plt.subplots(2, 2,figsize=(12,6))

for data in range(0,len(temp_list)): 
    whisker_count=0 
    if(data%2==0 and data!=0):
        whisker_count+=1
    i = []
    s = []
    r = []
    l = []  
    ax1 = axes[whisker_count][0]
    ax2=axes[whisker_count][1]
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

    for v in range(0,len(r)):
        total_force.append((-r[v]+l[v])/2)

    strain_value = []
    force_strain=[]
    for j in range(100,len(s),200):
        strain_value.append(s[j])
        force_strain.append(total_force[j])
   
    
    ax1.plot(i,total_force)
    ax2.plot(strain_value,force_strain, '-o')
    ax1.legend([f'Small whisker','Large whisker'])
    ax2.legend([f'Small whisker','Large whisker'])  
    
    ax2.set_xlabel('Strain')
    ax1.set_xlabel('Timestep (fs)')
    ax1.set_ylabel('Force (eV/Å)')
    ax2.set_ylabel('Force (eV/Å)')
    ax2.set_title(f'Force vs strain at {temp_list[data]} K')
    ax1.set_title(f'Force vs timestep at {temp_list[data]} K')
    ax1.grid()
    ax2.grid()
path = "plot_code/milestone_plots/"
save_path = os.path.join(path,f"m9_temp0_plot_small.png")
plt.savefig(save_path)
plt.tight_layout()
plt.show()