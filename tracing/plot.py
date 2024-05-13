import matplotlib as mat
import matplotlib.pyplot as plt
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

# Get Data
sphereSamples = readDataPoint("results/sphereSamples.csv")
# hits = readData("hits.csv")
# hits2 = [int(x)/sampleCount for x in hits["hits"]]
# orgs = [sphereSamples[int(i)] for i in hits["org_id"]]

hitsConst = readData("results/hitsConstant.csv")
hits2Const = [int(x)/sampleCount for x in hitsConst["hits"]]
orgs = [sphereSamples[int(i)] for i in hitsConst["org_id"]] # Should be the same

# Plotting
# Plot shadowing function
if(sphereSampler == "simple" or sphereSampler == "simple_quad"):
	orghits_by_phi = {}
	for i,org in enumerate(orgs):
		if org.phi not in orghits_by_phi:
			orghits_by_phi[org.phi] = ([],[], [])
		orghits_by_phi[org.phi][0].append(math.pi/2 - org.lat) # change to altitude
		# orghits_by_phi[org.phi][1].append(hits2[i])
		orghits_by_phi[org.phi][2].append(hits2Const[i])

	for phi,orghits in orghits_by_phi.items():
		# plt.plot(orghits[0], orghits[1], "-") # variable height
		plt.plot(orghits[0], orghits[2], "--", label="phi="+str(phi)) # constant height
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
