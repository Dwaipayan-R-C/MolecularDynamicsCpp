import matplotlib.pyplot as plt
import csv


i = []
s = []
r = []
l = []


with open('../data/milestone9.dat', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    next(plots)
    for row in plots:
        i.append(float(row[0]))
        s.append(float(row[1]))
        r.append(float(row[2]))
        l.append(float(row[3]))

plt.figure(0)

plt.title()
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('milestone_plots/m9plot.png')
i