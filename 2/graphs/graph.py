import os
os.chdir("D:\\СЛАУ\\2")

import matplotlib.pyplot as plt

f = open("data.txt")
a = f.readlines()
n = len(a)
d = []
t1 = []
t2 = []
for i in range (n):
    s = a[i].split()
    d.append(int(s[0]))
    t1.append(float(s[1]) * 1000)
    t2.append(float(s[2]) * 1000)
    
plt.figure(figsize=(16, 9))
plt.plot(d, t1, '-o', color = "blue", label = "full_matrix")
plt.plot(d, t2, '-o', color = "green", label = "CSR_matrix")
plt.legend()
plt.xlabel('N, количество ненулевых элементов')
plt.ylabel('T, мc')
plt.xlim(0, 8.1 * 10 ** 6)
plt.ylim(0, 11.5)
plt.grid(True, which='major', linewidth=1)
plt.grid(True, which='minor', linewidth=0.3)
plt.minorticks_on()
plt.savefig('graph.jpg', dpi = 300)
plt.show()