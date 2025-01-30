import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

def readCSVProbs(name, probIndex):
    probs = []
    reflectIndex = []
    reflectIndexStr = []
    with open(name, "r") as file:
        lines = file.readlines()
        
        # start = 0
        # for i,l in enumerate(lines):
        #     if(l.strip()==""):
        #         start = i + 2 ## header strip
        #         break
        
        # lines = lines[start:]
        
        found = False
        end = 0
        for i,l in enumerate(lines):
            if(l.strip()==""):
                end = i
                found = True
                break
        if found:
            lines = lines[0:end]

        for line in lines[1:]:
            column = line.split(",")
            reflectIndex.append(int(column[0]))
            reflectIndexStr.append(column[1])
            probs.append(float(column[probIndex]))
    return (reflectIndex, reflectIndexStr, probs)

def readCSVRebounceProbs(name, probIndex):
    probs = []
    with open(name, "r") as file:
        lines = file.readlines()

        found = False
        end = 0
        for i,l in enumerate(lines):
            if(l.strip()==""):
                end = i
                found = True
                break
        if found:
            lines = lines[0:end]

        for line in lines[1:]:
            column = line.split(",")
            probs.append(float(column[probIndex]))
    return probs

if len(sys.argv) != 3:
    raise ValueError("Need 2 runtime arguments: identifier sampleCount")
identifier = sys.argv[1]
sampleCount = sys.argv[2]

reflectIndexP, reflectIndexStrP, distP = readCSVProbs("./temp/dist_"+identifier+"_P.csv", 3)
reflectIndexL, reflectIndexStrL, distL = readCSVProbs("./temp/dist_"+identifier+"_L.csv", 3)
reflectIndexSpecular, reflectIndexStrSpecular, distSpecular = readCSVProbs("./temp/dist_"+identifier+"_S.csv", 3)
specularRebounceProb = readCSVRebounceProbs("./temp/dist_"+identifier+"_S.csv", 4)

plt.figure(figsize=(19.20,9.83))
ax = plt.subplot()
plt.bar(reflectIndexP, distP, 0.25, align="edge", label="BRDF_Paul", zorder=2)
# plt.bar(reflectIndexL, distL, 0.25, align="edge", label="BRDF_Lyanne", zorder=2)
plt.bar(reflectIndexSpecular, distSpecular, -0.25, align="edge", label="Mesh"+"("+sampleCount+"samples)")

plt.bar(reflectIndexSpecular, specularRebounceProb, -0.25, align="edge", bottom= np.array(distSpecular)-np.array(specularRebounceProb), label="rebounce")

plt.xticks(reflectIndexSpecular, rotation=90)
ax.set_xticklabels(reflectIndexStrSpecular)
plt.legend()
plt.title("Exit Probability Distribution ("+identifier+")")
plt.xlabel("Reflect Path")
plt.ylabel("Probability")
# plt.show()

fileName = "./images/dist_"+identifier+"_S"+str(sampleCount)+".png"
filePath = Path(fileName)
filePath.parent.mkdir(exist_ok=True, parents=False) # Ensure folder exists

plt.savefig(fileName, bbox_inches="tight", dpi=150)