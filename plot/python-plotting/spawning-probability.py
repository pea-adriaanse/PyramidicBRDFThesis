import matplotlib.pyplot as plt
import numpy as np

rho = 0.6
alpha = 54.7 * np.pi / 180

xEnd = 10

xs = np.linspace(0, xEnd,100)
ys = 1-np.exp(-rho*xs)

plt.xlim([0,xEnd])
plt.ylim([0,1])
plt.plot(xs, ys)
# plt.title("Spawning Probability")
plt.xlabel("$A$")
plt.ylabel("$P_C(A)$")
plt.grid()
plt.show()