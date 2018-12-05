#####
#   Linear Least Square Approximation
#####

from warnings import warn
import matplotlib.pyplot as plt
import numpy as np
from sympy import *
from Ploy_divide import extended_synthetic_division


def __expression():
    x, y, z = symbols("x y z")
    exprsn = sympify(input("Enter an expression : "))
    ev = float(input("Enter: "))
    f = lambdify(x, exprsn, 'numpy')
    print(exprsn)
    print(f(ev))


def __read_matrix(m, n):
    line = input("Matrix: ")
    a = []
    for i in line.replace(" ", "").split(","):
        a.append(float(i))
    return np.array(a).reshape(m, n), len(a)


def __sign(x):
    r = np.real(x)
    i = np.imag(x)
    if abs(x) == 0:
        return "+"
    elif r !=0 and i!=0:
        return "+"
    elif r >0 and i==0:
        return "+"
    elif r>0 and i>0:
        return "+"
    elif r==0 and i>0:
        return "+"
    elif r<0 and i>0:
        return

    elif np.real(x) >0:
        return "+"
    elif np.imag(x) > 0:
        return "+"
    else:
        return "-"


def __construct_poly(Coefficients):
    order = len(Coefficients)
    sign = ""

    equation = ""
    for i in range(order):
        equation += "+({})*x^{}".format(Coefficients[i], order - i - 1)
    return equation


if __name__ == "__main__":
    div, n, rem = extended_synthetic_division([1, 0, 0, 0, -6], [1, -3j])
    print(div)
    print(__construct_poly(div))
