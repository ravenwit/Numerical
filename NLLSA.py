#####
#   Non-Linear Least Square Approximation
#####

import matplotlib.pyplot as plt
import numpy as np
import math
import random
from scipy.linalg.blas import dgemm


def translate(value, leftMin, leftMax, rightMin, rightMax):
    # Figure out how 'wide' each range is
    if value > leftMax:
        value = leftMax
    elif value < leftMin:
        value = leftMin
    leftSpan = leftMax - leftMin
    rightSpan = rightMax - rightMin

    # Convert the left range into a 0-1 range (float)
    valueScaled = float(value - leftMin) / float(leftSpan)

    # Convert the 0-1 range into a value in the right range.
    return rightMin + (valueScaled * rightSpan)


def _solve_matrix(a, b, residuals=False):
    i = dgemm(alpha=1.0, a=a.T, b=a.T, trans_b=True)
    x = np.linalg.solve(i, dgemm(alpha=1.0, a=a.T, b=b)).flatten()

    if residuals:
        return x, np.linalg.norm(np.dot(a, x) - b)
    else:
        return x


"""
Functions 
"""


def func_gauss(x):
    return 10 * math.exp(-x ** 2 / (2 * 5 ** 2))


def func_exp(x):
    return math.exp(x)


def func_sin(x):
    return math.sin(x)


def func_line(x):
    return 5 * x + 8


def __generate_data_points(func, deviation):
    x_min = -100
    x_max = 100
    points = 400
    step = (x_max - x_min) / points
    x = np.full(points, 0.0)
    x[0] = x_min
    x = np.arange(x_min, x_max, step)
    y = np.array(list(map(lambda x: (func(x) + (random.random() * deviation) % deviation), x)))
    return x, y


def get_data():
    x = np.array([1, 2, 3, 4])
    y = np.array([6, 5, 7, 10])
    return x, y


def get_coeffs(x, y, o):
    Vandermonde = np.array(np.ones(len(x)))
    for i in range(1, o):
        __x = np.array(list(map(lambda x: pow(x, i), x)))
        Vandermonde = np.vstack([Vandermonde, __x])
    Vandermonde = Vandermonde.T
    wA = np.vstack([list(map(lambda t: (Vandermonde[t] * y[t]), range(len(y))))])
    wy = y * y.T
    # f_coeff  = _solve_matrix(wA, wy)
    f_coeff = _solve_matrix(Vandermonde, y)
    return f_coeff


def plot_data(x, coeff, step):
    order = len(coeff)
    x_map = np.arange(x[0], x[len(x) - 1], step)
    y_map = np.full(len(x_map), coeff[0])
    for i in range(len(x_map)):
        for j in range(1, order):
            y_map[i] += coeff[j] * x_map[i] ** j
    print(coeff)
    point_size = translate(len(y), 20, 100, 16, 1)
    plt.plot(x, y, 'o', label='Original data', markersize=point_size)
    line = plt.plot(x_map, y_map, 'r', label='Fitted line')
    plt.legend()
    plt.show()
    return line


if __name__ == "__main__":
    x, y = __generate_data_points(func_sin, 0.00001)
    while 1:
        order = int(input("Order of the curve: ")) + 1
        smootheness = 100
        plot_data(x, get_coeffs(x, y, order), 1 / smootheness)
