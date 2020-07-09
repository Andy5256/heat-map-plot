#coding=utf-8
import matplotlib.pyplot as plt
import numpy as np
points = np.arange(-5,5,0.01)
xs,ys = np.meshgrid(points,points)
z = np.sqrt(xs**2 + ys**2)
fig = plt.figure()
ax = fig.add_subplot(221)
ax.imshow(z)
ax = fig.add_subplot(222)
ax.imshow(z,cmap=plt.cm.gray)
ax = fig.add_subplot(223)
ax.imshow(z,cmap=plt.cm.cool)
ax = fig.add_subplot(224)
ax.imshow(z,cmap=plt.cm.hot)
plt.show()
