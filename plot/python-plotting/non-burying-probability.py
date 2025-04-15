import matplotlib.pyplot as plt
import numpy as np

rho = 0.6
alpha = 54.7 * np.pi / 180

xs = np.linspace(0, 2.5,100)
ys = np.exp(-rho * ( (2 * -xs)/(np.tan(alpha)) )**2)

plt.xlim([0,2.5])
plt.ylim([0,1])
plt.plot(xs, ys)
# plt.title("Non-Burying Probability")
plt.xlabel("$z$")
plt.ylabel("$P_E(z)$")
plt.grid()
plt.show()