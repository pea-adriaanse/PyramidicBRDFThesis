import matplotlib as mat
import matplotlib.pyplot as plt
import math
import csv

# Settings
with open("settings.csv", "r") as file:
	settings = {}
	for line in csv.reader(file):
		settings[line[0]] = line[1]

sphereSampleCount = int(settings["sphereSampleCount"])
peakSampleCount = int(settings["peakSampleCount"])

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

def get_org(org_id, sphereSampleCount):
	from collections import namedtuple
	Point = namedtuple('Point', 'x y z phi lat')
	longitude = org_id * math.pi/(1+math.sqrt(5.0))
	latitude = math.asin(org_id* 2/(2*sphereSampleCount+1))
	x = math.cos(longitude) * math.cos(latitude)
	y = math.sin(longitude) * math.cos(latitude)
	z = math.sin(latitude)
	return Point(x, y, z, longitude, latitude)

hits = readData("hits.csv")
hits2 = [int(x)/peakSampleCount for x in hits["hits"]]
# print(hits)

lats = [get_org(int(x), sphereSampleCount).lat for x in hits["org_id"]]
plt.plot(lats, hits2)
plt.show()