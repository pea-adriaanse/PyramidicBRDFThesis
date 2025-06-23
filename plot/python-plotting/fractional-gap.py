import matplotlib.pyplot as plt
import numpy as np

rho = 0.6
alpha = 54.7 * np.pi / 180

xEnd = 6

xs = np.linspace(0, xEnd, 200)
gap = np.exp(-4*rho*xs**2/np.tan(alpha)**2)
print(np.exp(-4*rho*10**2/np.tan(alpha)**2))

plt.xlim([0,xEnd])
plt.ylim([0,1])
plt.plot(xs, gap, label="$A_G'(h)$")
# plt.title("Fractional Gap")
plt.xlabel("$h\ [\mu m]$")
plt.legend()
plt.grid()
plt.savefig("fractional-gap-notitle.png", bbox_inches='tight')
plt.show()

ys = np.linspace(0, 1, 1000)
h = np.tan(alpha)*np.sqrt(-np.log(ys)*4*rho)
print(np.tan(alpha)*np.sqrt(-np.log(0.0000000000001)*4*rho))
plt.xlim([0,1])
plt.ylim([0,xEnd])
plt.plot(ys, h, label="$h_G(A_G')$")
plt.xlabel("$A_G'$")
plt.legend()
plt.grid()
plt.savefig("fractional-gap-height-notitle.png", bbox_inches='tight')
plt.show()