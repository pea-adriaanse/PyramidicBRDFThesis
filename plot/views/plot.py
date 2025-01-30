import matplotlib.pyplot as plt

ax = plt.figure().add_subplot(projection="3d")
ax.plot_surface(x,y,z)
ax.set(xlim=(xmin, xmax), ylim=(ymin,ymax), rcount=len(x), ccount=len(x[0]))
plt.show()
