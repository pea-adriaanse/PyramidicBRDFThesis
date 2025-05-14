import matplotlib.pyplot as plt
import numpy as np

rho = 0.6
alpha = 54.7 * np.pi / 180

xEnd = 6

xs = np.linspace(0, xEnd, 1000)
gap = np.exp(-4*rho*xs**2/np.tan(alpha)**2)

plt.xlim([0,xEnd])
plt.ylim([0,1])
plt.plot(xs, gap, label="$A_G'(h)$")
# plt.title("Fractional Gap")
plt.xlabel("$h$")
plt.legend()
plt.grid()
plt.show()

ys = np.linspace(0, 1, 1000)
h = np.tan(alpha)*np.sqrt(-np.log(ys)*4*rho)
plt.xlim([0,1])
plt.ylim([0,xEnd])
plt.plot(ys, h, label="$h_G(A_G')$")
plt.xlabel("$A_G'$")
plt.legend()
plt.grid()
plt.show()