import matplotlib as mat
import matplotlib.pyplot as plt
import numpy as np
import math
import csv

# Settings
with open("results/settings.csv", "r") as file:
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

# Plotting

# Plot landscape and highpoint & misses
landData = readData("results/land.csv")
fig = plt.figure()
ax = plt.axes(projection="3d")
landX = np.array([float(x) for x in landData["x"]])
landY = np.array([float(x) for x in landData["y"]])
landZ = np.array([float(x) for x in landData["z"]])
land = [x.reshape(landWidthSampleCount, landWidthSampleCount) for x in (landX, landY, landZ)]

cmap = plt.cm.plasma
norm=plt.Normalize(vmin=min(landZ),vmax=max(landZ))
ax.plot_surface(land[0], land[1], land[2],cmap=cmap, linewidth=0.1, zorder=1)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("Land Sampling")

ax.set_xlim(-width/2, width/2)
ax.set_ylim(-width/2, width/2)
ax.set_zlim(min(landZ), maxConstHeight+width)

heightPoints = readData("results/highPointSampling.csv")
# for i in range(len(heightPoints["orgx"])):
# 	(x,y,z) = tuple([float(heightPoints["org"+x][i]), float(heightPoints["hit"+x][i])] for x in ("x","y","z"))
# 	color = "red" if bool(heightPoints["hit"][i]) else "blue"
# 	ax.plot(x,y,z, c=color, zorder=4)
plt.show()

# Plot lyannes shadowing
shadow_lyanne = readData("results/shadow_lyanne.csv")
shadow_lyanne_theta = [float(x) for x in shadow_lyanne["theta"]]
shadow_lyanne_shadow = [float(x) for x in shadow_lyanne["shadow"]]
plt.plot(shadow_lyanne_theta, shadow_lyanne_shadow, "r-", label="lyanne")

# Plot shadowing function
sphereSamples = readDataPoint("results/sphereSamples.csv")
# hits = [int(x)/sampleCount for x in hitsData["hits"]]
# orgs = [sphereSamples[int(i)] for i in hitsData["org_id"]]
hitsConstData = readData("results/hitsConstant.csv")
hitsConst = [int(x)/sampleCount for x in hitsConstData["hits"]]
orgs = [sphereSamples[int(i)] for i in hitsConstData["org_id"]] # Should be the same

heights = [int(x) for x in readData("results/hitsHeightConstant.csv")["occurrance"]]

colormap = plt.cm.cool
norm=plt.Normalize(vmin=0,vmax=math.pi/2)
if(sphereSampler == "simple" or sphereSampler == "simple_quad"):
	orghits_by_phi = {}
	for i,org in enumerate(orgs):
		if org.phi not in orghits_by_phi:
			orghits_by_phi[org.phi] = ([],[], [])
		orghits_by_phi[org.phi][0].append(math.pi/2 - org.lat) # change to altitude
		# orghits_by_phi[org.phi][1].append(hits[i])
		orghits_by_phi[org.phi][2].append(hitsConst[i])

	for phi,orghits in orghits_by_phi.items():
		# plt.plot(orghits[0], orghits[1], "-") # variable height
		plt.plot(orghits[0], orghits[2], "--",color=colormap(norm(phi)), label="phi="+str(phi)) # constant height
	plt.legend()
	# plt.legend(["variable height", "constant height"])

else:
	lats = [x.lat for x in orgs]
	# latsConst = [x.lat for x in orgsConst]
	plt.plot(lats, hits2, label="variable heights")
	plt.plot(lats, hits2Const, label="constant heights")

plt.axvline((90-54.7) * math.pi / 180)
plt.xscale("linear")
plt.title("Shadowing Function")
plt.xlabel("latitude")
plt.ylabel("visibility")
plt.show()

# Plot shadowing function by height
hitsConstByHeightData = readData("results/hitsConstantByHeight.csv")
hitsConstByHeight = [float(x) for x in hitsConstByHeightData["hits"]]
hitsConstByHeightPolar = [math.pi/2-float(x) for x in hitsConstByHeightData["lat"]]
length = len(hitsConstByHeightData["hits"])
thetaCount = int(length/heightBins)

colormap = plt.cm.cool
norm=plt.Normalize(vmin=0,vmax=heightBins)
for heightBin in range(heightBins):
	start = heightBin*thetaCount
	stop = (heightBin+1)*thetaCount
	x = hitsConstByHeightPolar[start:stop]
	y = hitsConstByHeight[start:stop]
	plt.plot(x, y, color=colormap(norm(heightBin)), label="heightBin="+str(heightBin))

plt.legend()
plt.axvline((90-54.7) * math.pi / 180)
plt.xscale("linear")
plt.title("Shadowing Function by Height")
plt.xlabel("latitude")
plt.ylabel("visibility")
plt.show()

# Plot hits by height
heights = readData("results/hitsHeightConstant.csv")
heightsBins = len(heights["occurrance"])
plt.plot([int(x) for x in heights["occurrance"]], "-*", label="occurrance")
plt.plot([int(x) for x in heights["hit"]], "--*", label="hits")
plt.xlim(0, heightsBins)
plt.xlabel("bin")
plt.ylabel("count")
plt.legend()
plt.title("Height Dependance")
plt.show()

# Plot hits by height ratio
heightsRelative = [int(heights["hit"][i])/int(heights["occurrance"][i]) if int(heights["occurrance"][i]) != 0 else None for i in range(heightsBins)]
plt.plot(heightsRelative, "-*")
plt.xlabel("bin")
plt.ylabel("hit ratio")
plt.title("Height Dependance")
plt.show()

# Plot height distribution
sampleHeights = readData("results/heightDistributionConstant.csv");
plt.title("Height Sample Distribution")
plt.ylabel("count")
plt.xlabel("bin")
plt.plot([float(x) for x in sampleHeights["count"]], "-*")
plt.show()

# Plot sphere samples
fig = plt.figure()
ax = plt.axes(projection="3d")
ax.scatter([x.x for x in orgs], [x.y for x in orgs], [x.z for x in orgs])
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("Sphere Sampling")
plt.show()
