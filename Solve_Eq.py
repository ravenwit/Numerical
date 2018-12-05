import numpy as np
import math
import functools
import cmath
from scipy.optimize import *
import sympy as sp
from sympy import I


def __expression():
    x, y, z = sp.symbols("x y z")
    exprsn = sp.sympify(input("Enter a function : "))
    f = sp.lambdify(x, exprsn, 'numpy')
    f = functools.partial(f)
    return exprsn, f


def __extended_synthetic_division(dividend, divisor):
    '''Fast polynomial division by using Extended Synthetic Division. Also works with non-monic polynomials.'''
    # dividend and divisor are both polynomials, which are here simply lists of coefficients. Eg: x^2 + 3x + 5 will be represented as [1, 3, 5]

    out = list(dividend)  # Copy the dividend
    normalizer = divisor[0]
    for i in range(len(dividend) - (len(divisor) - 1)):
        out[i] /= normalizer  # for general polynomial division (when polynomials are non-monic),
        # we need to normalize by dividing the coefficient with the divisor's first coefficient
        coef = out[i]
        if coef != 0:  # useless to multiply if coef is 0
            for j in range(1,
                           len(divisor)):  # in synthetic division, we always skip the first coefficient of the divisor,
                # because it is only used to normalize the dividend coefficients
                out[i + j] += -divisor[j] * coef

    # The resulting out contains both the quotient and the remainder, the remainder being the size of the divisor (the remainder
    # has necessarily the same degree as the divisor since it is what we couldn't divide from the dividend), so we compute the index
    # where this separation is, and return the quotient and remainder.
    separator = -(len(divisor) - 1)
    return out[:separator], out[separator:]  # return quotient, remainder.


def __construct_poly(Coefficients):
    order = len(Coefficients)
    equation = ""
    for i in range(order):
        equation += "+({})*x^{}".format(Coefficients[i], order - i - 1)
    return equation


def __func(f, x):
    res = f(x)
    return res


def __func_d(f, x):
    h = 0.000001
    res = (f(x + h) - f(x)) / h
    return res


def __func_d_d(f, x):
    h = 0.000001
    res = (f(x + h) - 2 * f(x) + f(x - h)) / h ** 2
    return res


def bisection(f):
    lim_a = 0
    lim_b = 0
    while True:
        lim_a = float(input("Enter lower limit: "))
        lim_b = float(input("Enter upper limit: "))
        if __func(f, lim_a) * __func(f, lim_b) < 0: break
    std = 0.001
    while True:
        c = (lim_b + lim_a) / 2
        if __func(f, c) * __func(f, lim_a) < 0:
            lim_b = c
        elif __func(f, c) * __func(f, lim_b) < 0:
            lim_a = c
        if abs(lim_b - lim_a) <= std:
            print(c)
            print(__func(f, c))
            break


def newton(f, xk):
    std = 0.0000000000000000001
    while (abs(__func(f, xk))) >= std:
        xk -= __func(f, xk) / __func_d(f, xk)
    print(xk, __func(f, xk))


def laguerre(f, n, xk):
    std = 1e-5
    while (abs(__func(f, xk))) >= std:
        G = __func_d(f, xk) / __func(f, xk)
        H = G ** 2 - __func_d_d(f, xk) / __func(f, xk)
        dnmtr_p = G + cmath.sqrt((n - 1) * (n * H - G ** 2))
        dnmtr_n = G - cmath.sqrt((n - 1) * (n * H - G ** 2))
        dnmtr = dnmtr_p if abs(dnmtr_p) > abs(dnmtr_n) else dnmtr_n
        a = n / dnmtr
        xk -= a
        # if abs(xk) - np.real(xk) ==0: xk = np.real(xk)
    # print(xk, __func(f, xk))
    print(xk)
    return xk


def find_all_root(expression):
    x = sp.Symbol("x")
    f_poly = sp.Poly(expression)
    f_poly = f_poly.subs(I, 1j) if len(f_poly.monoms()[0]) >1 else f_poly
    coeffs = f_poly.all_coeffs()
    order = f_poly.degree()
    f = sp.lambdify(x, expression, 'numpy')
    roots = []
    roots.append(laguerre(f, order, (5+1j) if len(roots) == 0 else roots[-1]))
    div_poly_coeffs, rem = __extended_synthetic_division(coeffs, [1, (-1)*roots[-1]])
    if len(div_poly_coeffs) > 1:
        expression = __construct_poly(div_poly_coeffs)
        find_all_root(expression)
    return roots


if __name__ == "__main__":
    # while True:
    #     sfunc, f = __expression()
    #     x = int(input("Enter point: "))
    #     print(__func_d(f, x))
    sfunc, f = __expression()
    # n = int(input("Enter order: "))
    # laguerre(f, n, 3)
    # newton(f, 0.22252094585341875+0.974927927510545j)
    find_all_root(sfunc)
