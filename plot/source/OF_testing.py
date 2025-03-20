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
    return R(L,T) + A(T) * math.exp(- pow( ((L-L0(T))**2)/(2*(S(T)**2)), 1.5))


xs = np.arange(300, 1000, 1)

for t in [0,40,60,80]:
    y = [magic(x, t) for x in xs]
    plt.plot(xs, y)
plt.show()