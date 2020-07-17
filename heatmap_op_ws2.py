#coding=utf-8
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.transforms as mtransforms
points = np.arange(-5,5,0.1)
xs,ys = np.meshgrid(points,points)
# z1 = np.loadtxt("shiftvec_333_100_op_hdch.dat")
# z = np.loadtxt("shiftvec_111_100_op_4315.dat") # sc for GeS
z = np.loadtxt("shiftvec_ws2_222_op.dat") # sc for WS2
# z = np.transpose(z1)+z2
nkp = 100
N = np.zeros((2*nkp,2*nkp))
for i in range(nkp) :
    for j in range(nkp) :
        N[i,j] = N[i+nkp,j] = N[i,j+nkp] = N[i+nkp,j+nkp] = z[i,j]
# N = np.transpose(N)  # nkps in order of gamma Y ; X M
fig = plt.figure()
ax = fig.add_subplot(111)
# ax.imshow(z,cmap=plt.cm.seismic,interpolation='bicubic',vmin=-100,vmax=100)
def do_plot(ax, Z, transform):
    im = ax.imshow(Z,cmap=plt.cm.summer,interpolation='none',vmin=-10,vmax=10)
    trans_data = transform + ax.transData
    im.set_transform(trans_data)

    # display intended extent of the image
    # x1, x2, y1, y2 = im.get_extent()
    # ax.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], "y--",
    #         transform=trans_data)
    ax.set_xlim(-5, 150)
    ax.set_ylim(-5, 150)
# image skew
do_plot(ax, z, mtransforms.Affine2D().skew_deg(15, 15))
# # ax = fig.add_subplot(222)
# # ax.imshow(z,cmap=plt.cm.seismic,vmin=-10,vmax=10)
# # ax = fig.add_subplot(223)
# # ax.imshow(z,cmap=plt.cm.seismic,vmin=-10,vmax=10)
# # ax = fig.add_subplot(224)
# # ax.imshow(z,cmap=plt.cm.seismic,vmin=-10,vmax=10)
plt.show()
# plt.savefig("sv_333_100_op_hdch.png")
# plt.close()