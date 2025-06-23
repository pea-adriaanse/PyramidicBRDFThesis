import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import sys

arguments = sys.argv
assert len(arguments) == 4

dir1 = arguments[1]
dir2 = arguments[2]
destination = arguments[3]

def readArrayFile(file):
	with open(file, 'r') as f:
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
# thetas_x = readArrayFile("thetas.txt")
# phis_y = readArrayFile("phis.txt")
# errors_z = np.array(readArrayFile("errors.txt"))
# exitProbBRDF_z = np.array(readArrayFile("exitProbBRDF.txt"))
# exitProbMesh_z = np.array(readArrayFile("exitProbMesh.txt"))

resolution = 64

thetas_x_1 = (readArrayFile(dir1+"thetas.txt"))
thetas_x_2 = (readArrayFile(dir2+"thetas.txt"))
assert thetas_x_1 == thetas_x_2
phis_x_1 = (readArrayFile(dir1+"phis.txt"))
phis_x_2 = (readArrayFile(dir2+"phis.txt"))
assert phis_x_1 == phis_x_2
errors_z_1 = np.array(readArrayFile(dir1+"errors.txt"))
errors_z_2 = np.array(readArrayFile(dir2+"errors.txt"))
errors_z_diff = errors_z_1 - errors_z_2

X, Y = np.meshgrid(thetas_x_1, phis_x_1)
ax = plt.figure().add_subplot(projection='3d')
ax.plot_surface(X,Y,errors_z_diff, rcount=resolution, ccount=resolution, cmap=mpl.colormaps["plasma"])
ax.set_xlabel("theta")
ax.set_ylabel("phi")
# ax.set(xlim=xBounds[0:2],ylim=yBounds[0:2])
ax.set_zlim([-0.1, 0.6])
# ax.set_zlim([0, 0.15])
ax.view_init(azim=ax.azim-75)
plt.savefig(destination)
plt.show()

errors_z_1_flat = errors_z_1.flatten()
errors_z_2_flat = errors_z_2.flatten()

print(errors_z_1_flat.sum())
print(errors_z_2_flat.sum())

# plt.hist(errors_z_1_flat, weights=np.ones_like(errors_z_1_flat)/len(errors_z_1_flat), bins=50, range=(0,1), histtype='step', label="paul")
# plt.hist(errors_z_2_flat, weights=np.ones_like(errors_z_2_flat)/len(errors_z_2_flat), bins=50, range=(0,1), histtype='step', label="lyanne")
plt.hist(errors_z_1_flat, bins=40, range=(0,1), histtype='step', label="paul")
plt.hist(errors_z_2_flat, bins=40, range=(0,1), histtype='step', label="lyanne")
plt.show()

# fig = plt.figure()
# # ax = fig.add_subplot(1,2,1,projection='3d')
# ax = fig.add_subplot(projection='3d')
# min_val = 0.9 #exitProbBRDF_z.min()
# max_val = 1 #(1-min_val)*1.2+min_val
# ax.set_zlim([min_val,max_val])
# ax.plot_surface(X,Y,exitProbBRDF_z, rcount=resolution, ccount=resolution, cmap=mpl.colormaps["plasma"], vmin=min_val, vmax= max_val)
# ax.set_xlabel("theta")
# ax.set_ylabel("phi")
# ax.view_init(azim=ax.azim-75)
# plt.savefig(dir+"totalExitProbBRDF.png")
# plt.show()

# # ax = fig.add_subplot(1,2,2,projection='3d')
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# # use same min_val & max_val. (may fail)
# # ax.set_zlim([0.7,1])
# ax.set_zlim([min_val,max_val])
# ax.plot_surface(X,Y,exitProbMesh_z, rcount=resolution, ccount=resolution, cmap=mpl.colormaps["plasma"], vmin=min_val, vmax= max_val)
# ax.set_xlabel("theta")
# ax.set_ylabel("phi")
# ax.view_init(azim=ax.azim-75)
# plt.savefig(dir+"totalExitProbMesh.png")
# plt.show()