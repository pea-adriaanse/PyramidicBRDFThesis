import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

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
thetas_x = readArrayFile("./thetas.txt")
phis_y = readArrayFile("./phis.txt")
errors_z = np.array(readArrayFile("./errors.txt"))
exitProbBRDF_z = np.array(readArrayFile("./exitProbBRDF.txt"))
exitProbMesh_z = np.array(readArrayFile("./exitProbMesh.txt"))

X, Y = np.meshgrid(thetas_x, phis_y)
ax = plt.figure().add_subplot(projection='3d')
ax.plot_surface(X,Y,errors_z, rcount=64, ccount=64, cmap=mpl.colormaps["plasma"])
ax.set_xlabel("theta")
ax.set_ylabel("phi")
# ax.set(xlim=xBounds[0:2],ylim=yBounds[0:2])
plt.show()

ax = plt.figure().add_subplot(projection='3d')
min_val = exitProbBRDF_z.min()
max_val = (1-min_val)*1.2+min_val
ax.set_zlim([0.7,1])
ax.plot_surface(X,Y,exitProbBRDF_z, rcount=64, ccount=64, cmap=mpl.colormaps["plasma"], vmin=min_val, vmax= max_val)
ax.set_xlabel("theta")
ax.set_ylabel("phi")
plt.show()

ax = plt.figure().add_subplot(projection='3d')
# use same min_val & max_val. (may fail)
ax.set_zlim([0.7,1])
ax.plot_surface(X,Y,exitProbMesh_z, rcount=64, ccount=64, cmap=mpl.colormaps["plasma"], vmin=min_val, vmax= max_val)
ax.set_xlabel("theta")
ax.set_ylabel("phi")
plt.show()