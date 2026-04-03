import os
import math
os.chdir("D:\\eque\\HW\\graphs")

import matplotlib.pyplot as plt

f = open("iter_data_sym.txt")
a = f.readlines()
n = len(a)
d = []
r1 = []
r2 = []
r3 = []
for i in range (n):
    v1, v2, v3 = list(map(float, a[i].split()))
    if v1 > 0:
        r1.append(math.log(v1))
    if v2 > 0:
        r2.append(math.log(v2))
    if v3 > 0:
        r3.append(math.log(v3))
    d.append(i + 1)
    
plt.figure(figsize=(16, 9))
plt.plot(d[:len(r1)], r1, '-o', markersize=3, color = "blue", label = "Jacobi_speed")
plt.plot(d[:len(r2)], r2, '-o', markersize=3, color = "green", label = "Gauss_speed")
plt.plot(d[:len(r3)], r3, '-o', markersize=3, color = "purple", label = "Gauss_sym_speed")
plt.legend()
plt.title("Зависимость невязки от итерации")
plt.xlabel('n, итерация')
plt.ylabel('ln(Невязка)')
plt.grid(True, which='major', linewidth=1)
plt.grid(True, which='minor', linewidth=0.3)
plt.minorticks_on()
plt.savefig('graph_iterations.jpg', dpi = 300)
plt.show()

plt.figure(figsize=(16, 9))

f = open("time_jac_speed.txt")
a = f.readlines()
n = len(a)
t1 = [0] * n
r1 = [0] * n
for i in range (n):
    t1[i], r1[i] = list(map(float, a[i].split()))
for i in range (n):
    r1[i] = math.log(abs(r1[i]))
plt.plot(t1, r1, '-o', markersize=3, color = "blue", label = "Jacobi_speed")

f = open("time_gauss_speed.txt")
a = f.readlines()
n = len(a)
t2 = [0] * n
r2 = [0] * n
for i in range (n):
    t2[i], r2[i] = list(map(float, a[i].split()))
for i in range (n):
    r2[i] = math.log(abs(r2[i]))
plt.plot(t2, r2, '-o', markersize=3, color = "green", label = "Gauss_sym_speed")

f = open("time_gauss_sym_speed.txt")
a = f.readlines()
n = len(a)
t3 = [0] * n
r3 = [0] * n
for i in range (n):
    t3[i], r3[i] = list(map(float, a[i].split()))
for i in range (n):
    r3[i] = math.log(abs(r3[i]))
plt.plot(t3, r3, '-o', markersize=3, color = "Yellow", label = "Gauss")

plt.legend()
plt.title("Зависимость невязки от времени")
plt.xlabel('t, мс')
plt.ylabel('ln(Невязка)')
plt.grid(True, which='major', linewidth=1)
plt.grid(True, which='minor', linewidth=0.3)
plt.minorticks_on()
plt.savefig('graph_time.jpg', dpi = 300)
plt.show()