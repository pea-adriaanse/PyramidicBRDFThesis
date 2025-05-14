import matplotlib.pyplot as plt
import numpy as np

rho = 0.6
alpha = 54.7 * np.pi / 180

xEnd = 2.5

xs = np.linspace(0, xEnd,100)
cdf = 1-np.exp(-rho*4*xs**2/np.tan(alpha)**2)
pdf = 8*rho*xs/(np.tan(alpha)**2) * np.exp(-rho*4*xs**2/np.tan(alpha)**2)

plt.xlim([0,xEnd])
plt.ylim([0,1])
plt.plot(xs, cdf, label="$P_B(z)$")
plt.plot(xs, pdf, label="$E(z)$")
# plt.title("Spawning Probability")
plt.xlabel("$z$")
plt.legend()
plt.grid()
plt.show()