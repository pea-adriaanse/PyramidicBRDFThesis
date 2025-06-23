import matplotlib.pyplot as plt
import matplotlib as mpl

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
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
plot = ax.plot_surface(xs,ys,zs, rcount=len(xs), edgecolor='0.7', linewidth=0.1)#, cmap=mpl.colormaps["plasma"])
# plot = ax.plot_surface(xs,ys,zs, rcount=len(xs), edgecolor='0.7', linewidth=0.1, cmap=mpl.colormaps["plasma"])
# fig.colorbar(plot, shrink=0.8, pad=0.1)
ax.set(xlim=xBounds[0:2],ylim=yBounds[0:2])
ax.set_xlabel("$\\theta$")
ax.set_ylabel("$\phi$")
print(type(zLimitsAutomatic))
print(zLimitsAutomatic)
if not zLimitsAutomatic:
	ax.set_zlim(zLimits)
# ax.view_init(azim=ax.azim-75)
# print(zs.max())
# ax.set_zlim(0, 1)
if len(savePlotPath)>0:
	plt.savefig(savePlotPath,bbox_inches='tight')
plt.show()

# plt.figure()
# # plt.contourf(xs, ys, zs)
# plt.pcolor(xs,ys,zs, cmap=mpl.colormaps["plasma"])
# plt.colorbar()
# plt.xlabel("$theta$")
# plt.ylabel("$phi$")
# plt.show()

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