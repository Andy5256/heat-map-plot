#coding=utf-8
import matplotlib.pyplot as plt
import numpy as np
points = np.arange(-2,2,0.1)
xs,ys = np.meshgrid(points,points)
# print(xs)
z = np.loadtxt("shiftvec_333_40.dat")
print(np.size(xs,1))
print(np.size(z,1))
fig = plt.figure()
ax = fig.add_subplot(111)
# plt.plot(z)
ax.imshow(z,cmap=plt.cm.seismic)
# plt.show()
plt.savefig("sv_333_wt_40.png")
plt.close()