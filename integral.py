#!/usr/bin/env python

import math
import time


def __func(x):
    # return x*math.exp(x**2)
    # return (x*math.exp(math.asin(x)))/float(math.sqrt(1-x**2))
    return (x ** 4) * math.exp(-x)
    # return math.sin(x)
    # return x**2


def __func_d(x):
    # return (2*x**2+1)*math.exp(x**2)
    # return x**2*math.exp(math.asin(x))/(math.sqrt(1-x**2)*(1-x**2))
    # +x*math.exp(math.asin(x))/(1-x**2)+math.exp(math.asin(x))/math.sqrt(1-x**2)
    return (4 * x ** 3) * math.exp(-x) - (x ** 4) * math.exp(-x)
    # return math.cos(x)
    # return 2*x


def __func_k(x, x_k):
    return (x - x_k) * __func_d(x_k) + __func(x_k)


def _integral_interpolate(lim_a, lim_b, step):
    sum = 0
    x_k = lim_a
    x_k1 = lim_a + step
    s_t = time.clock()
    while x_k1 <= lim_b:
        sum += step * (__func_k(x_k, x_k) + __func_k(x_k1, x_k))
        x_k += step
        x_k1 += step
    e_t = time.clock()
    return sum / 2, (e_t - s_t)


def _integral_trapezoid(lim_a, lim_b, step):
    sum = 0
    i = lim_a
    s_t = time.clock()
    while i + step <= lim_b:
        sum += __func(i) + __func(i + step)
        i += step
    e_t = time.clock()
    return (step / 2) * sum, (e_t - s_t)


def _integral_simpson(lim_a, lim_b, step):
    sum = 0
    i = lim_a / 2
    s_t = time.clock()
    while (2 * (i + step)) <= lim_b:
        sum += __func(2 * i) + 4 * __func(2 * i + step) + __func(2 * (i + step))
        i += step
    e_t = time.clock()
    return (step / 3) * sum, (e_t - s_t)


if __name__ == "__main__":
    func = "x^2"
    print("\nThis program serves three method for integration \n (1) Trapezoid Method\n "
          "(2) Simpson Method\n (3) Interpolation Method\n")
    method = int(input("\n\tWhich method do you wanna use to integrate {} (1/2/3): ".format(func)))
    limit_l = float(input("\tEnter lower limit: "))
    limit_h = float(input("\tEnter upper limit: "))
    part = int(input("\tEnter division number: "))
    _h = (limit_h - limit_l) / part
    print("\n\tStep : {}".format(_h))
    if method is 1:
        res = _integral_trapezoid(limit_l, limit_h, _h)
    elif method is 2:
        res = _integral_simpson(limit_l, limit_h, _h)
    else:
        res = _integral_interpolate(limit_l, limit_h, _h)
    print("\n\tResult : {}".format(res[0]))
    print("\n\tTime needed : {} seconds".format(res[1]))
