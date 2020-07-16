#coding=utf-8
import matplotlib.pyplot as plt
import numpy as np
points = np.arange(-5,5,0.1)
xs,ys = np.meshgrid(points,points)
# z1 = np.loadtxt("shiftvec_333_100_op_hdch.dat")
z = np.loadtxt("shiftvec_111_100_op_4315.dat")
# z = np.transpose(z1)+z2
nkp = 100
N = np.zeros((2*nkp,2*nkp))
for i in range(nkp) :
    for j in range(nkp) :
        N[i,j] = N[i+nkp,j] = N[i,j+nkp] = N[i+nkp,j+nkp] = z[i,j]
# N = np.transpose(N)  # nkps in order of gamma Y ; X M
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(N,cmap=plt.cm.seismic,vmin=-100,vmax=100)
# # ax = fig.add_subplot(222)
# # ax.imshow(z,cmap=plt.cm.seismic,vmin=-10,vmax=10)
# # ax = fig.add_subplot(223)
# # ax.imshow(z,cmap=plt.cm.seismic,vmin=-10,vmax=10)
# # ax = fig.add_subplot(224)
# # ax.imshow(z,cmap=plt.cm.seismic,vmin=-10,vmax=10)
plt.show()
# plt.savefig("sv_333_100_op_hdch.png")
# plt.close()