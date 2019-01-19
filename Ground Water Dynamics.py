import numpy as np
import matplotlib as mp
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter



nx = 100
ny = 50
ni = 5

sigma0 = 1
a = -0.04
phi0 = 200
b = -20
lx = 1000
hx = lx / nx
ly = 500
hy = ly / ny

phi = np.array(np.zeros((nx + 1, ny + 1)))
sigma = np.array(np.zeros((nx + 1, ny + 1)))
f = np.array(np.zeros((nx + 1, ny + 1)))

p = 1.2


def relaxation(p, hx, hy, u, d, s):
    h2 = hx * hx
    a = h2 / (hy * hy)
    b = 1 / (4 * (1 + a))
    ab = a * b
    q = 1 - p
    for i in range(1, nx):
        for j in range(1, ny):
            xp = b * (d[i + 1][j] / d[i][j] + 1)
            xm = b * (d[i - 1][j] / d[i][j] + 1)
            yp = ab * (d[i][j + 1] / d[i][j] + 1)
            ym = ab * (d[i][j - 1] / d[i][j] + 1)
            u[i][j] = q * u[i][j] + p * (xp * u[i + 1][j]
                                         + xm * u[i - 1][j] + yp * u[i][j + 1]
                                         + ym * u[i][j - 1] + h2 * s[i][j])
            return u


for i in range(0, nx + 1):
    x = i * hx
    for j in range(0, ny + 1):
        y = j * hy
        sigma[i][j] = sigma0 + a * ny
        phi[i][j] = phi0 + b * np.cos(np.pi * x / lx) * y / ly
        f[i][j] = 0

for step in range(0, ni):
    for j in range(0, ny):
        phi[0][j] = (4 * phi[1][j] - phi[2][j]) / 3
        phi[nx][j] = (4 * phi[nx - 1][j] - phi[nx - 2][j]) / 3

    phi = relaxation(p, hx, hy, phi, sigma, f)


fig = plt.figure()
ax = fig.gca(projection='3d')

X = np.arange(0, ly+hy, hy)
Y = np.arange(0, lx+hx, hx)
X, Y = np.meshgrid(X, Y)

surf= ax.plot_surface(X, Y, phi, cmap=cm.coolwarm,
                linewidth=0, antialiased=False)


ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)


# fig.show()
plt.show()
