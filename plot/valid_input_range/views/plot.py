import matplotlib.pyplot as plt

# FINAL CONTOOUR:
# ax = plt.figure().add_subplot()
# contours = ax.contourf(xs,ys,zs, levels=[-0.5,1], colors=["red"])
# artists, labels = contours.legend_elements()
# plt.xlabel("$\\theta$ [rad]")
# plt.ylabel("$\\phi$ [rad]")
# plt.title("Backbounce input domain")
# plt.legend(artists, ["domain"])
# ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
# plt.show()

# PLOT 3D
ax = plt.figure().add_subplot(projection='3d')
ax.plot_surface(xs,ys,zs, rcount=len(xs), )
ax.set(xlim=xBounds[0:2],ylim=yBounds[0:2])
ax.set_xlabel("$theta$")
ax.set_ylabel("$phi$")
# ax.set_zlim(-1, 1)
plt.show()

# ax = plt.figure().add_subplot(projection='3d')
# ax.plot_surface(xs,ys,zs, rcount=len(xs)/100, ccount=len(xs[0])/100)
# ax.set(xlim=xBounds, ylim=yBounds)
# ax.set_xlabel("X")
# ax.set_ylabel("Y")
# plt.show()

# import numpy as np
# delta = 0.025
# x = np.arange(-3.0, 3.0, delta)
# y = np.arange(-2.0, 2.0, delta)
# X, Y = np.meshgrid(x, y)
# Z1 = np.exp(-X**2 - Y**2)
# Z2 = np.exp(-(X - 1)**2 - (Y - 1)**2)
# Z = (Z1 - Z2) * 2
# print(Z)

# print(">>>>>>>")

# print(zs)

# ax.contour(X, Y, Z)
# plt.show()