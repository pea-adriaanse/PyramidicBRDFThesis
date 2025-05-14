import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mat
import matplotlib.pyplot as plt
import numpy as np
import math
import csv

rho = 0.6
alpha = 54.7 * np.pi / 180

def DE(theta, phi):
	return np.sin(alpha)*np.sin(theta)*np.cos(phi)+np.cos(alpha)*np.cos(theta)

def DW(theta, phi):
	return -np.sin(alpha)*np.sin(theta)*np.cos(phi)+np.cos(alpha)*np.cos(theta)

def DN(theta, phi):
	return np.sin(alpha)*np.sin(theta)*np.sin(phi)+np.cos(alpha)*np.cos(theta)

def DS(theta, phi):
	return -np.sin(alpha)*np.sin(theta)*np.sin(phi)+np.cos(alpha)*np.cos(theta)

def Dot(theta, phi):
	return max(0,DE(theta,phi))+max(0,DW(theta,phi))+max(0,DN(theta,phi))+max(0,DS(theta,phi))

def G1(theta,phi):
	return 4.0*np.cos(alpha)*np.cos(theta)/Dot(theta,phi)

def GL(theta):
	if(theta <= np.pi/2-alpha):
		return 1
	return 4*np.cos(alpha)*np.cos(theta)/(2*np.cos(alpha)*np.cos(theta)+np.cos(theta-alpha))

xEnd = np.pi/2

xs = np.linspace(0, xEnd, 100)
gLs = [GL(x) for x in xs]
# plt.plot(xs, gLs, 'r--', dashes=(3,3), linewidth=3, label="$G_L$", zorder=3)

plt.xlim([0,xEnd])
plt.ylim([0,1.2])

# for phi in np.linspace(0, np.pi/4, 4):
# 	phi = round(phi, 2)
# 	g1s = [G1(x, phi) for x in xs]
# 	plt.plot(xs, g1s, '-.', label="$G_1(\phi="+str(phi)+")$")

g1s_0 = [G1(x, 0) for x in xs]
plt.plot(xs, g1s_0, 'r-', label="$G_1(\phi=0)$", zorder=2)

g1s_pi_4 = [G1(x, np.pi/4) for x in xs]
plt.plot(xs, g1s_pi_4, 'b-', label="$G_1(\phi=\pi/4)$", zorder=1)
plt.vlines([np.pi/2-alpha], [0],[1.2])

# IMPORTED OLD PYTHON PLOTTING:

# Settings
with open("tracing/results/settings.csv", "r") as file:
	settings = {}
	for line in csv.reader(file):
		settings[line[0]] = line[1]

sphereSampleCount = int(settings["sphereSampleCount"])
sampleCount = int(settings["sampleCount"])
sphereSampler = settings["sphereSampler"]
landWidthSampleCount = int(settings["landWidthSampleCount"])
width = float(settings["width"])
maxConstHeight = float(settings["maxConstHeight"])
heightBins = int(settings["heightBins"])
simpleSphereLongCount = int(settings["simpleSphereLongCount"])
simpleSphereLatCount = int(settings["simpleSphereLatCount"])

# Utility
def readData(fileName):
	with open(fileName, "r") as file:
		first = file.readline().split('=')
		sep = ','
		if first[0] == 'sep':
			sep = first[1]
		else:
			file.seek(0)

		header = next(file).strip().split(sep)
		data = {}
		for h in header:
			data[h] = []

		for line in csv.reader(file, delimiter=sep):
			for i, l in enumerate(line):
				data[header[i]].append(l)
	return data

def readDataPoint(fileName):
	with open(fileName, "r") as file:
		first = file.readline().split('=')
		sep = ','
		if first[0] == 'sep':
			sep = first[1]
		else:
			file.seek(0)

		header = next(file).strip().split(sep)
		from collections import namedtuple
		Point = namedtuple('Point', header)

		data = []

		for line in csv.reader(file, delimiter=sep):
			lineNum = [float(x) for x in line]
			data.append(Point(*lineNum))
	return data

# Plot shadowing function
sphereSamples = readDataPoint("tracing/results/sphereSamples.csv")
# hits = [int(x)/sampleCount for x in hitsData["hits"]]
# orgs = [sphereSamples[int(i)] for i in hitsData["org_id"]]
hitsConstData = readData("tracing/results/hitsConstant.csv")
hitsConst = [int(x)/sampleCount for x in hitsConstData["hits"]]
hitsVarData = readData("tracing/results/hitsVar.csv")
hitsVar = [int(x)/sampleCount for x in hitsVarData["hits"]]
orgs = [sphereSamples[int(i)] for i in hitsConstData["org_id"]] # Should be the same

heights = [int(x) for x in readData("tracing/results/hitsHeightConstant.csv")["occurrance"]]
heightsVar = [int(x) for x in readData("tracing/results/hitsHeightVar.csv")["occurrance"]]

colormap = plt.cm.cool
norm=plt.Normalize(vmin=0,vmax=math.pi/2)
orghits_by_phi = {}
orghits_by_phiVar = {}
for i,org in enumerate(orgs):
	if org.phi not in orghits_by_phi:
		orghits_by_phi[org.phi] = ([],[], [])
	orghits_by_phi[org.phi][0].append(math.pi/2 - org.lat) # change to altitude
	orghits_by_phi[org.phi][1].append(hitsVar[i])
	orghits_by_phi[org.phi][2].append(hitsConst[i])

i = 0
for phi,orghits in orghits_by_phi.items():
	if(i >= 2):
		break # Plot only 0 and pi/4
	i = i + 1
	if(i==1):
		phiName = "0"
		color = "r"
	else:
		phiName="\pi/4"
		color = "b"
	# print([x > 1 for x in orghits[2]])
	# print([x > 1 for x in orghits[1]])
	# plt.plot(orghits[0], orghits[1], "-", label="variable height $\phi$="+phiName, zorder=2, linewidth=1) # variable height
	plt.plot(orghits[0], orghits[2], color+"--", label="$G_{MC}(\phi="+phiName+")$", zorder=1, dashes=(4,4), linewidth=3) # constant height
	# plt.plot(orghits[0], orghits[2], "--",color=colormap(norm(phi)), label="phi="+str(phi)) # constant height
# plt.legend(["variable height", "constant height"])
plt.legend()



# plt.title("Spawning Probability")
plt.ylabel("visibility")
plt.xlabel("$\\theta$")
plt.legend()
plt.grid()
plt.show()