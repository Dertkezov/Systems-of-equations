import os
import math
os.chdir("D:\\eque\\HW\\graphs")

import matplotlib.pyplot as plt

f = open("iter_data.txt")
a = f.readlines()
n = len(a)
d = []
r1 = []
r2 = []
r3 = []
r4 = []
r5 = []
r6 = []
for i in range (n):
    v1, v2, v3, v4, v5, v6 = list(map(float, a[i].split()))
    if v1 > 0:
        r1.append(math.log(v1))
    if v2 > 0:
        r2.append(math.log(v2))
    if v3 > 0:
        r3.append(math.log(v3))
    if v4 > 0:
        r4.append(math.log(v4))
    if v5 > 0:
        r5.append(math.log(v5))
    if v6 > 0:
        r6.append(math.log(v6))
    d.append(i + 1)
    
plt.figure(figsize=(16, 9))
plt.plot(d[:len(r1)], r1, '-o', markersize=3, color = "blue", label = "Gauss")
plt.plot(d[:len(r2)], r2, '-o', markersize=3, color = "green", label = "SOR")
plt.plot(d[:len(r3)], r3, '-o', markersize=3, color = "yellow", label = "Gauss with acceleration")
plt.plot(d[:len(r4)], r4, '-o', markersize=3, color = "purple", label = "Fast_slope")
plt.plot(d[:len(r5)], r5, '-o', markersize=3, color = "red", label = "GC")
plt.plot(d[:len(r6)], r6, '-o', markersize=3, color = "magenta", label = "GMRES")
plt.legend()
plt.title("Зависимость невязки от итерации")
plt.xlabel('n, итерация')
plt.ylabel('ln(Невязка)')
plt.grid(True, which='major', linewidth=1)
plt.grid(True, which='minor', linewidth=0.3)
plt.minorticks_on()
plt.savefig('graph_iterations_ellips_matrix.jpg', dpi = 300)
plt.show()

plt.figure(figsize=(16, 9))

f = open("time_gauss.txt")
a = f.readlines()
n = len(a)
t1 = [0] * n
r1 = [0] * n
for i in range (n):
    t1[i], r1[i] = list(map(float, a[i].split()))
for i in range (n):
    r1[i] = math.log(abs(r1[i]))
plt.plot(t1, r1, '-o', markersize=3, color = "blue", label = "Gauss")

f = open("time_sor.txt")
a = f.readlines()
n = len(a)
t2 = [0] * n
r2 = [0] * n
for i in range (n):
    t2[i], r2[i] = list(map(float, a[i].split()))
for i in range (n):
    r2[i] = math.log(abs(r2[i]))
plt.plot(t2, r2, '-o', markersize=3, color = "green", label = "SOR")

f = open("time_chebyshev.txt")
a = f.readlines()
n = len(a)
t3 = [0] * n
r3 = [0] * n
for i in range (n):
    t3[i], r3[i] = list(map(float, a[i].split()))
for i in range (n):
    r3[i] = math.log(abs(r3[i]))
plt.plot(t3, r3, '-o', markersize=3, color = "Yellow", label = "Gauss_with_acceleration")

f = open("time_slope.txt")
a = f.readlines()
n = len(a)
t4 = [0] * n
r4 = [0] * n
for i in range (n):
    t4[i], r4[i] = list(map(float, a[i].split()))
for i in range (n):
    r4[i] = math.log(abs(r4[i]))
plt.plot(t4, r4, '-o', markersize=3, color = "purple", label = "Fast_descent")

f = open("time_cg.txt")
a = f.readlines()
n = len(a)
t5 = [0] * n
r5 = [0] * n
for i in range (n):
    t5[i], r5[i] = list(map(float, a[i].split()))
for i in range (n):
    r5[i] = math.log(abs(r5[i]))
plt.plot(t5, r5, '-o', markersize=3, color = "red", label = "GC")

f = open("time_gm.txt")
a = f.readlines()
n = len(a)
t6 = [0] * n
r6 = [0] * n
for i in range (n):
    t6[i], r6[i] = list(map(float, a[i].split()))
for i in range (n):
    r6[i] = math.log(abs(r6[i]))
plt.plot(t6, r6, '-o', markersize=3, color = "magenta", label = "GMRES")

plt.legend()
plt.title("Зависимость невязки от времени")
plt.xlabel('t, мс')
plt.ylabel('ln(Невязка)')
plt.grid(True, which='major', linewidth=1)
plt.grid(True, which='minor', linewidth=0.3)
plt.minorticks_on()
plt.savefig('graph_time_ellips_matrix_1.jpg', dpi = 300)
plt.show()

plt.figure(figsize=(16, 9))
plt.plot(t3, r3, '-o', markersize=3, color = "Yellow", label = "Gauss_with_acceleration")
plt.plot(t4, r4, '-o', markersize=3, color = "purple", label = "Fast_descent")
plt.plot(t5, r5, '-o', markersize=3, color = "red", label = "GC")
plt.plot(t6, r6, '-o', markersize=3, color = "magenta", label = "GMRES")
plt.legend()
plt.title("Зависимость невязки от времени")
plt.xlabel('t, мс')
plt.ylabel('ln(Невязка)')
plt.grid(True, which='major', linewidth=1)
plt.grid(True, which='minor', linewidth=0.3)
plt.minorticks_on()
plt.savefig('graph_time_ellips_matrix_2.jpg', dpi = 300)
plt.show()