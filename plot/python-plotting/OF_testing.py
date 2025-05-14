import matplotlib.pyplot as plt
import numpy as np
import math

def lerp(x, L, R):
    return x*(R-L)+L

def R(L,T):
    if T < 40:
        return 0.06
    return lerp(T/40-1,0.06, 0.5)

def A(T):
    if T < 40:
        return 0.84
    return lerp(T/40-1, 0.84, 0.12)

def L0(T):
    if T < 40:
        return lerp(T/40, 590, 550)
    return lerp(T/40-1, 550, 500)

def S(T):
    return lerp(T/80, 40, 30)

def magic(L, T):
    T = min(80,T)
    return R(L,T) + A(T) * math.exp(- pow( ((L-L0(T))**2)/(2*(S(T)**2)), 1.5))

xlow = 380
xhigh = 780
xs = np.arange(xlow, xhigh, 1)

# ts = np.arange(0,90,1)
ts = [0,20,40,60,80]
colors = plt.colormaps["plasma"](np.linspace(0.1,0.9,len(ts)))

for i,t in enumerate(reversed(ts)):
    y = [magic(x, t) for x in xs]
    plt.plot(xs, y, label="$\phi="+str(t)+"\degree$", color=colors[i])
plt.xlim([xlow, xhigh])
plt.ylim([0,1])
plt.grid()
plt.legend()
plt.ylabel("reflectance")
plt.xlabel("wavelength [nm]")
plt.show()