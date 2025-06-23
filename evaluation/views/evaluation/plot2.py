import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import sys

arguments = sys.argv
dir = "./results/"
if len(arguments) == 2:
	dir = arguments[1]

def readArrayFile(file):
	with open(dir+file, 'r') as f:
		text = f.readline().strip().replace(' ', '') # Remove all whitespace on ends and spaces in between.
	return readArrayString(text)

def readArrayString(string):
	assert len(string)>0
	assert string[0] == '[' and string[-1]==']'
	if(string[1] == '['): # array of arrays
		arr = []
		while True:
			assert string[1] == '[' # all entries are arrays
			j = string.index(']')
			arr.append(readArrayString(string[1:j+1])) # read entry
			if string[j+1] == ',': # continues
				string = string[j+1:] # trim off last read entry (except ',')
			else:
				assert string[j+1] == ']' # end of array
				return arr
	else: # array of values
		return [float(x) for x in string[1:-1].split(',')]

# Context:
thetas_x = readArrayFile("thetas.txt")
phis_y = readArrayFile("phis.txt")
errors_z = np.array(readArrayFile("errors.txt"))
exitProbBRDF_z = np.array(readArrayFile("exitProbBRDF.txt"))
exitProbMesh_z = np.array(readArrayFile("exitProbMesh.txt"))

resolution = 64

X, Y = np.meshgrid(thetas_x, phis_y)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
min_val = 0#errors_z.min()
# max_val = errors_z.max()
max_val = float.__ceil__(errors_z.max()*10)/10
surface = ax.plot_surface(X,Y,errors_z, rcount=resolution, ccount=resolution, cmap=mpl.colormaps["plasma"], vmin=min_val, vmax=max_val)
ax.set_xlabel("$\\theta$")
ax.set_ylabel("$\phi$")
fig.colorbar(surface, ax=ax, shrink=0.8)
# ax.set(xlim=xBounds[0:2],ylim=yBounds[0:2])
# ax.set_zlim([0, 0.9])
ax.view_init(azim=ax.azim-75)
plt.savefig(dir+"absoluteErrorSum.png",bbox_inches='tight')
print(errors_z.max())
plt.show()

plt.hist(errors_z.flatten(), bins=50, range=(0,1))
plt.savefig(dir+"absoluteErrorSumHistogram.png",bbox_inches='tight')
plt.show()

fig = plt.figure()
# ax = fig.add_subplot(1,2,1,projection='3d')
ax = fig.add_subplot(projection='3d')
min_val = 0.9 #exitProbBRDF_z.min()
max_val = 1 #(1-min_val)*1.2+min_val
ax.set_zlim([min_val,max_val])
surface = ax.plot_surface(X,Y,exitProbBRDF_z, rcount=resolution, ccount=resolution, cmap=mpl.colormaps["plasma"], vmin=min_val, vmax= max_val)
fig.colorbar(surface, ax=ax, shrink=0.8)
ax.set_xlabel("$\\theta$")
ax.set_ylabel("$\phi$")
ax.view_init(azim=ax.azim-75)
plt.savefig(dir+"totalExitProbBRDF.png",bbox_inches='tight')
plt.show()

# ax = fig.add_subplot(1,2,2,projection='3d')
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
# use same min_val & max_val. (may fail)
# ax.set_zlim([0.7,1])
ax.set_zlim([min_val,max_val])
surface = ax.plot_surface(X,Y,exitProbMesh_z, rcount=resolution, ccount=resolution, cmap=mpl.colormaps["plasma"], vmin=min_val, vmax= max_val)
fig.colorbar(surface, ax=ax, shrink=0.8)
ax.set_xlabel("$\\theta$")
ax.set_ylabel("$\phi$")
ax.view_init(azim=ax.azim-75)
plt.savefig(dir+"totalExitProbMesh.png",bbox_inches='tight')
plt.show()