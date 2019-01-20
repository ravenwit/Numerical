#!/usr/bin/env python

"""
    Importing necessary libraries
"""
import numpy as np
import matplotlib as mp
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

###########################################################################################

"""
    Declaring necessary physical constants
    
"""

sigma0 = 1
a = -0.04
phi0 = 200
b = -20

# Length and width of boundary
lx = 1000
ly = 500

###########################################################################################


"""
    Setting up grid mesh

"""

nx = 1000
ny = 500
ni = 5

hx = lx / nx
hy = ly / ny

X = np.arange(0, ly + hy, hy)  # In case of Mesh Y and L_x are correlated
Y = np.arange(0, lx + hx, hx)  # In case of Mesh X and L_y are correlated
X, Y = np.meshgrid(X, Y)

# Optimising constant
omega = 2 / (1 + np.pi / np.sqrt(nx * ny))
# omega = 1.2
###########################################################################################


"""
    Initiating the HEAD FUNCTION, HYDRAULIC CONDUCTIVITY and INFILTRATION
"""

phi = np.array(np.zeros((nx + 1, ny + 1)))
sigma = np.array(np.zeros((nx + 1, ny + 1)))
f = np.array(np.zeros((nx + 1, ny + 1)))

for i in range(0, nx + 1):
    x = i * hx
    for j in range(0, ny + 1):
        y = j * hy
        sigma[i][j] = sigma0 + a * ny
        phi[i][j] = phi0 + b * np.cos(np.pi * x / lx) * y / ly
        f[i][j] = 0


###########################################################################################


def relaxation(func, conduct, infilt):
    '''
        Relaxation scheme with difference equation

    :param func: Function with corresponding boundary condition and arbitrary values for other points in spaces
    :param conduct: Conductivity in points of space
    :param infilt: Infiltration rate in points of space
    :return: Revised solution function
    '''

    global omega
    global nx
    global ny
    global hx
    global hy

    alpha = hx ** 2 / hy ** 2
    factor0 = 1 / (4 * (1 + alpha))
    factor1 = alpha * factor0
    p = 1 - omega

    for i in range(1, nx):
        for j in range(1, ny):
            term1 = factor0 * (conduct[i + 1][j] / conduct[i][j] + 1)
            term2 = factor0 * (conduct[i - 1][j] / conduct[i][j] + 1)
            term3 = factor1 * (conduct[i][j + 1] / conduct[i][j] + 1)
            term4 = factor1 * (conduct[i][j - 1] / conduct[i][j] + 1)

            #   Evaluating the function with the updated function as the previous function
            #   phi(k+1 th guess) = (1-p) * phi(kth guess) + p * phi(updated with difference equation)

            func[i][j] = p * func[i][j] + omega * (term1 * func[i + 1][j]
                                                   + term2 * func[i - 1][j] + term3 * func[i][j + 1]
                                                   + term4 * func[i][j - 1] + (hx ** 2) * infilt[i][j])
            return func


def plot_surface(X, Y, func):
    '''
        Plotting a surface curve with wireframe diagram
    :param X:   Mesh Grid
    :param Y:   Mesh Grid
    :param func:    2D array
    :return:
    '''

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.set_xlabel("X", fontsize=10)
    ax.set_ylabel("Y", fontsize=10)
    ax.set_zlabel("$\phi (x, y)$", fontsize=12, rotation=-90)

    ax.plot_wireframe(X, Y, func, label='Water Head')

    # plt.title('Water head surface')
    ax.legend()
    fig.show()


def plot_color_surface(X, Y, func):
    '''
        Plotting a surface with high low color scheme
    :param X:   Mesh Grid
    :param Y:   Mesh Grid
    :param func:    2D array
    :return:
    '''

    fig, ax = plt.subplots()
    ax = fig.gca(projection='3d')

    surf = ax.plot_surface(X, Y, phi, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    ax.set_xlabel("X", fontsize=10)
    ax.set_ylabel("Y", fontsize=10)
    ax.set_zlabel("$\phi (x, y)$", fontsize=12, rotation=-90)

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)

    # plt.title('Water head surface')
    ax.legend()
    fig.show()


if __name__ == '__main__':

    # Initialising relaxation scheme satisfying boundary condition
    for step in range(0, ni):
        for j in range(0, ny):
            phi[0][j] = (4 * phi[1][j] - phi[2][j]) / 3
            phi[nx][j] = (4 * phi[nx - 1][j] - phi[nx - 2][j]) / 3

        phi = relaxation(phi, sigma, f)

    plot_surface(X, Y, phi)
    plot_color_surface(X, Y, phi)
    plt.show()
