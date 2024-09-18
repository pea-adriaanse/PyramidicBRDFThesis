import numpy as np
import matplotlib.pyplot as plt

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

reflectIndex, reflectIndexStr, dist = readCSVProbs("../iridescence/dist.csv", 2)
reflectIndexSpecular, _, distSpecular = readCSVProbs("../tracing/results/distSpecular.csv", 3)

ax = plt.subplot()
plt.bar(reflectIndex, dist, 0.35, align="edge", label="BRDF");
plt.bar(reflectIndexSpecular, distSpecular, -0.35, align="edge", label="Mesh");
plt.xticks(reflectIndex, rotation=90)
ax.set_xticklabels(reflectIndexStr)
plt.legend()
plt.title("Exit Probability Distribution");
plt.xlabel("Reflect Path")
plt.ylabel("Probability")
plt.show()