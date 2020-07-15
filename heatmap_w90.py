#coding=utf-8
import matplotlib.pyplot as plt
import numpy as np
points = np.arange(-10,10,0.05)
xs,ys = np.meshgrid(points,points)
z = np.loadtxt("shiftvec_333_200.dat")
N = np.zeros((400,400))
for i in range(200) :
    for j in range(200) :
        N[i,j] = N[i+200,j] = N[i,j+200] = N[i+200,j+200] = z[i,j]
# print(N)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(N,cmap=plt.cm.seismic,vmin=-10,vmax=10)
# # ax = fig.add_subplot(222)
# # ax.imshow(z,cmap=plt.cm.seismic,vmin=-10,vmax=10)
# # ax = fig.add_subplot(223)
# # ax.imshow(z,cmap=plt.cm.seismic,vmin=-10,vmax=10)
# # ax = fig.add_subplot(224)
# # ax.imshow(z,cmap=plt.cm.seismic,vmin=-10,vmax=10)
# plt.show()
plt.savefig("sv_333_200_w90.png")
plt.close()
