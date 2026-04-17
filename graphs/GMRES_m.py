import os
import math
os.chdir("D:\\eque\\HW\\graphs")

import matplotlib.pyplot as plt

f = open("iter_m.txt")
a = f.readlines()
n = len(a)
m = []
r = []
for i in range (n):
    v1, v2 = list(map(float, a[i].split()))
    m.append(int(v1))
    r.append(v2)

plt.figure(figsize=(16, 9))
plt.plot(m, r, '-o', markersize=3, color = "green")
plt.xlabel('Порядок m')
plt.ylabel('Общее количество итераций')
plt.grid(True, which='major', linewidth=1)
plt.grid(True, which='minor', linewidth=0.3)
plt.minorticks_on()
plt.savefig('graph_GMRES_iter.jpg', dpi = 300)
plt.show()


f = open("time_m.txt")
a = f.readlines()
n = len(a)
m = []
r = []
for i in range (n):
    v1, v2 = list(map(float, a[i].split()))
    m.append(int(v1))
    r.append(v2)

plt.figure(figsize=(16, 9))
plt.plot(m, r, '-o', markersize=3, color = "green")
plt.xlabel('Порядок m')
plt.ylabel('Общее время, мс')
plt.grid(True, which='major', linewidth=1)
plt.grid(True, which='minor', linewidth=0.3)
plt.minorticks_on()
plt.savefig('graph_GMRES_time.jpg', dpi = 300)
plt.show()

f = open("time_m.txt")