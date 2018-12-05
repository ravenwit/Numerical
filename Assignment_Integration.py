#!/usr/bin/env python

import math
import time


def __func(x, f):
    """
    Evaluating functions for different points.
    Since there are three functions to consider, the integer argument 'f' is used
    to identify the functions sinx, cosx and x^2

    :param x: The point for evaluation
    :param f: Enum for identifying function
    :return:  The value of function at point x
    """

    if f is 0:
        return math.sin(x)
    elif f is 1:
        return math.cos(x)
    else:
        return x ** 2


def __func_d(x, f):
    """
    Evaluating the differentials of functions for different points.
    Since there are three functions to consider, the integer argument 'f' is used
    to identify the derivative functions of sinx, cosx and x^2

    :param x: The point for evaluation
    :param f: Enum for identifying function
    :return:  The value of derivative at point x
    """

    if f is 0:
        return math.cos(x)
    elif f is 1:
        return -math.sin(x)
    else:
        return 2 * x


def __func_i(lim_a, lim_b, f):
    """
    Evaluating the analytical integrals of functions between the given limit: lim_a and lim_b.
    This function is used for the purpose of determining the accuracy of different methods
    of integration.
    The accuracy is measured with respect to the analytical result returned by this function.

    Since there are three functions to consider, the integer argument 'f' is used
    to identify the derivative functions of sinx, cosx and x^2

    :param lim_a: Lower limit
    :param lim_b: Upper limit
    :param f: Enum for identifying function
    :return:  The value of integration at point x
    """

    if f is 0:
        return -(math.cos(lim_b) - math.cos(lim_a))
    elif f is 1:
        return math.sin(lim_b) - math.sin(lim_a)
    else:
        return (1 / 3) * (lim_b ** 3 - lim_a ** 3)


def __func_k(x, x_k, f):
    """
    This function is defined for the purpose of integration by interpolation method

    :param x:       A point for evaluation
    :param x_k:     The next point of x
    :param f:       Enum for identifying function
    :return:        The value of (x - x_k) * (d/dx(f(x_k))  + f(x_k)
    """

    return (x - x_k) * __func_d(x_k, f) + __func(x_k, f)


def _integral_interpolate(lim_a, lim_b, step, f):
    """
    Integration by method of interpolation

    :param lim_a:   Lower limit
    :param lim_b:   Upper limit
    :param step:    The infinitesimal part of respective axis, dx or h
    :param f:       Enum for identifying function
    :return:        The value integration of function from lim_a to lim_b and the time needed for the process
    """
    sum = 0
    x_k = lim_a
    x_k1 = lim_a + step
    s_t = time.clock()  # The time before starting the process
    while x_k1 <= lim_b:
        # sum += step * (__func_k(x_k, x_k, f) + __func_k(x_k1, x_k, f))
        sum += step * (__func(x_k, f) + __func_k(x_k1, x_k, f))
        x_k += step
        x_k1 += step
    e_t = time.clock()  # The time at the end of the process
    return float("{0:.11f}".format(sum / 2)), float("{0:.3f}".format(e_t - s_t))


def _integral_trapezoid(lim_a, lim_b, step, f):
    """
    Integration by trapezoid rule
    S = (h/2)*Sum( f(x_i) + f(x_(i+1) )

    :param lim_a:   Lower limit
    :param lim_b:   Upper limit
    :param step:    The infinitesimal part of respective axis, dx or h
    :param f:       Enum for identifying function
    :return:        The value integration of function from lim_a to lim_b and the time needed for the process
    """
    sum = 0
    i = lim_a
    s_t = time.clock()  # The time before starting the process
    while i + step <= lim_b:
        sum += __func(i, f) + __func(i + step, f)
        i += step
    e_t = time.clock()  # The time before starting the process
    return float("{0:.11f}".format((step / 2) * sum)), float("{0:.3f}".format(e_t - s_t))


def _integral_simpson(lim_a, lim_b, step, f):
    """
    Integration by Simpson rule
    S = (h/3)*Sum( f(2x_i) + 4*f(x_(2i+1)) + f(x_(2i+2)) )

    :param lim_a:   Lower limit
    :param lim_b:   Upper limit
    :param step:    The infinitesimal part of respective axis, dx or h
    :param f:       Enum for identifying function
    :return:        The value integration of function from lim_a to lim_b and the time needed for the process
    """
    sum = 0
    i = lim_a / 2
    s_t = time.clock()  # The time before starting the process
    while (2 * (i + step)) <= lim_b:
        sum += __func(2 * i, f) + 4 * __func(2 * i + step, f) + __func(2 * (i + step), f)
        i += step
    e_t = time.clock()  # The time before starting the process
    return float("{0:.11f}".format((step / 3) * sum)), float("{0:.3f}".format(e_t - s_t))


if __name__ == "__main__":
    """
        Main Function
    """
    func = ["Sinx", "Cosx", "X^2"]  # A list containing the name of the function for displaying
    print("\nIntegration of sinx, cosx and x^2 using trapezoid and simpson method")

    # Asking user for the limits
    limit_l = float(input("\n\tEnter lower limit: "))
    limit_h = float(input("\tEnter upper limit: "))

    # Dividing the area under the function in the range of limits given by 1000000
    part = 1000000

    # Calculating and displaying the dx or h
    _h = (limit_h - limit_l) / part
    print("\tUsing h = {}".format(_h))

    print("\n{:^5}{:^17}{:^7}{:^13}{:^17}{:^7}{:^13}{:^17}{:^7}{:^13}".format("", "Trapezoid", "Time(s)", "Accuracy(%)",
                                                                              "Simpson", "Time(s)", "Accuracy(%)",
                                                                              "Interpolation", "Time(s)",
                                                                              "Accuracy(%)"))
    """
    Declaring some lists to contain the results of integration.
    First three lists are 3x2 in size containing three elements for three functions 
    sinx, cosx and x^2.
    Each element contains the integration result and time required.
    
    The last list holds the accuracy of the results obtained.
    """
    res_trap = []
    res_simp = []
    res_itrp = []
    accuracy = []

    # Index i is used as the substitute of Enum 'f'
    for i in range(3):
        res_trap.append(_integral_trapezoid(limit_l, limit_h, _h, i))
        res_simp.append(_integral_simpson(limit_l, limit_h, _h, i))
        res_itrp.append(_integral_interpolate(limit_l, limit_h, _h, i))

    # Calculating accuracy of the integration results
    for i in range(3):
        analytic_res = __func_i(limit_l, limit_h, i)
        accuracy.append([float("{0:.9f}".format(100 - abs((res_trap[i][0] - analytic_res) * (100 / analytic_res)))),
                         float(
                             "{0:.9f}".format(100 - abs((res_simp[i][0] - analytic_res) * (100 / analytic_res)))),
                         float(
                             "{0:.9f}".format(100 - abs((res_itrp[i][0] - analytic_res) * (100 / analytic_res))))])

    # Displaying the results in standard output
    for i in range(3):
        print("{:^5}{:^17}{:^7}{:^13}{:^17}{:^7}{:^13}{:^17}{:^7}{:^13}".format(func[i], res_trap[i][0], res_trap[i][1],
                                                                                accuracy[i][0], res_simp[i][0],
                                                                                res_simp[i][1], accuracy[i][1],
                                                                                res_itrp[i][0], res_itrp[i][1],
                                                                                accuracy[i][2]))
