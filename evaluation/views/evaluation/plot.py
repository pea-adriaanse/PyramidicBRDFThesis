import matplotlib.pyplot as plt
import numpy as np

# Context:
# thetas_x
# phis_y
# errors_z
# exitProbBRDF_z
# exitProbMesh_z

X, Y = np.meshgrid(thetas_x, phis_y)
ax = plt.figure().add_subplot(projection='3d')
ax.plot_surface(X,Y,errors_z) # rcount
ax.set_xlabel("theta")
ax.set_ylabel("phi")
# ax.set(xlim=xBounds[0:2],ylim=yBounds[0:2])
plt.show()

ax.plot_surface(X,Y,exitProbBRDF_z)
ax.set_xlabel("theta")
ax.set_ylabel("phi")
plt.show()

ax.plot_surface(X,Y,exitProbMesh_z)
ax.set_xlabel("theta")
ax.set_ylabel("phi")
plt.show()