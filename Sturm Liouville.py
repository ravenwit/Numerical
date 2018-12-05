#!/usr/bin/env python

"""
    Plot Legendre Equation with eigenvalues 0,1,2,3,4,5,6,7 in range -1,1 with conditions 1,1,0 (or) 0,0,1,1
    Plot Bessel Equation with eigenvalues 2 (or) 3 (or) 4 (or) 5 in range 0,80 with conditions 0.00000001,0,1

    Example:
        [1] Legendre Equation
	    [2] Bessel Equation
	    [3] Custom Equation

	    Select (1/2/3): 1

	    Enter some eigenvalues (with comma to separate): 0,1,2,3,4,5,6,7
	    Enter plot range (two values of x separated with comma): -1,1
	    Conditions (with comma to separate)
        Format: initial parameter, initial value, boundary parameter, boundary value. Example: 0,0,1,1
        Or initial parameter, initial function value, initial derivative value. Example: 1,1,0
        Enter conditions: 1,1,0

    Another Example:
        [1] Legendre Equation
	    [2] Bessel Equation
	    [3] Custom Equation

	    Select (1/2/3): 2

	    Enter some eigenvalues (with comma to separate): 3
	    Enter plot range (two values of x separated with comma): 0,150
	    Conditions (with comma to separate)
        Format: initial parameter, initial value, boundary parameter, boundary value. Example: 0,0,1,1
        Or initial parameter, initial function value, initial derivative value. Example: 1,1,0
        Enter conditions: 0.000001,0,1

"""

"""
    Importing necessary libraries
"""
import sys
import copy
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

###########################################################################################


"""
    Declaring some lists for Butcher tableau and F_i holds the weight coefficients for Runge-Kutta method
    All of these variables are necessary for Runge-Kutta method. 
"""
a = []
b = []
c = []
F_i = []
###########################################################################################

"""
    The functions for Sturm Liouville equation.
    They are globally declared for programming convenience. It's not required anymore to get
    them into functions argument.
"""
p = 0
p_d = 0
q = 0
s = 0
###########################################################################################

"""
#   Conditions is a list of the condition values of differential equation acquired from condition
    equation. User needs to find that out themselves and provide when program asks for it.
 
 
#   PROB is a string, being either 'IVP' (Initial value problem) or 'BVP' (Boundary value problem)
    depending upon the condition provided for the Sturm-Liouville Differential Equation.
 
    If the conditions list holds three value indicating the initial parameter and initial values 
    of the solution function and the 1st derivative of the solution function (Sturm-Liouville 
    equation is a second order differential equation, thus higher order of derivatives at the 
    initial value of parameter is not required for it to be a initial value problem as well as 
    solving for the function).


#   step is number to add or subtract from the current value of parameter to get the next or 
    previous value of the solution function respectively. 
    x + dx
    step = dx
"""
conditions = 0

PROB = ""

step = 0.01


###########################################################################################


def __set_conditions():
    """
    This function asks for the condition for solving the SL Diff Eq from the user.
    Once given it checks if it is a valid IVP or BVP condition.
    If not it exits the program.
    Otherwise it updates the global conditions list and global PROB storing proper
    values.

    :return: Nothing
    """
    global PROB
    global conditions
    conditions = list(map(float, [i for i in input(
        "Conditions (with comma to separate)\n"
        "Format: initial parameter, initial value, boundary parameter, boundary value. Example: 0,0,1,1\n"
        "Or initial parameter, initial function value, initial derivative value. Example: 1,1,0\n"
        "Enter conditions: ").split(',')]))
    cond_type = len(conditions)

    #   If user provides three number, then it's a IVP
    if cond_type == 3:
        PROB = 'IVP'
    #   If user provides four number, then it's a BVP
    elif cond_type == 4:
        PROB = 'BVP'
    else:
        print("Invalid conditions...")
        sys.exit()


def Butcher_tableau():
    """
        This function just inserts the value of Butcher Matrix to its
        elements a, b, c.


    :return: Returns the arrays a, b, c. Which isnt necessary. a, b and c are global
             variable.
    """

    # Referring to the global versions of a, b, c and F_i
    global a
    global b
    global c
    global F_i

    # If method equals 1 then the order is 4
    root21 = np.sqrt(21)
    c = np.array([0, 1, 1 / 2, 2 / 3, (7 - root21) / 14, (7 + root21) / 14, 1])
    a = np.array([
        [0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0],
        [3 / 8, 1 / 8, 0, 0, 0, 0, 0],
        [8 / 27, 2 / 27, 8 / 27, 0, 0, 0, 0],
        [3 * (3 * root21 - 7) / 392, 8 * (7 - root21) / 392, 48 * (7 - root21) / 392, 3 * (21 - root21) / 392, 0, 0, 0],
        [-5 * (231 + 51 * root21) / 1960, -40 * (7 + root21) / 1960, -320 * root21 / 1960,
         3 * (21 + 121 * root21) / 1960, 392 * (6 + root21) / 1960, 0, 0],
        [15 * (22 + 7 * root21) / 180, 120 / 180, 40 * (7 * root21 - 5) / 180, -63 * (3 * root21 - 2) / 180,
         -14 * (49 + 9 * root21) / 180, 70 * (7 - root21) / 180, 0]
    ])
    b = np.array([9 / 180, 0, 64 / 180, 0, 49 / 180, 49 / 180, 9 / 180])

    F_i = np.array(np.zeros(len(b)))

    return a, b, c


def __func_y1(y1, y2, eig_val, x):
    """
        This function gives the first derivative of the solution function.

        By theoretical definition y1 is the solution function and y2 is the
        first derivative of the solution function.

        The need for y1, x and eig_val in this function is that we will later use a general
        indicator for both __func_y1 and __func_y2. Therefore, their arguments must
        have the same form.

    :param y1: The solution function
    :param y2: First derivative of the solution function
    :param x:  Parameter
    :return: y2, First derivative of the solution function
    """
    return y2


def __func_y2(y1, y2, eig_val, x):
    """
        This function gives the second derivative of the solution function.

        By theoretical definition y1 is the solution function and y2 is the
        first derivative of the solution function.


    :param y1: The solution function
    :param y2: First derivative of the solution function
    :param x:  Parameter
    :return: Second derivative of the solution function
    """
    global p
    global p_d
    global q
    global s
    global conditions

    return s(x) - (1 / p(x)) * (p_d(x) * y2 + q(x, eig_val) * y1)


def F(func, y1, y2, eig_val, x):
    """
        Calculating the value of F_i for a single time step.

        We are actually solving two differential equations here to obtain a solution
        for the 2nd order differential.

    :param func: Which function to use, __func_y1 or __func_y2.
    :param y1: The solution function
    :param y2: First derivative of the solution function
    :param eig_val: Given eigenvalue
    :param x: Parameter
    :return: returns nothing. But update the global variable F_i for latter use.
    """
    # Referring to the global version of F_i
    global F_i
    # Calculating according to the theory. The article has a detailed explanation.
    for i in range(len(b)):
        F_i[i] = step * func((y1 + sum([a[i][j] * F_i[j] for j in range(len(a))])),
                             (y2 + sum([a[i][j] * F_i[j] for j in range(len(a))])),
                             eig_val, x + step * c[i])


def _runge_kutta(init_cond, end_boundary, eig_val, points, further):
    """
        The Runge-Kutta method to solve differential equation.
        This function can be regarded as the kernel to the whole algorithm.

    :param init_cond: Initial conditions
    :param end_boundary: The final boundary parameter. Necessary for solving BVP to calculate
                            only up-to that point.
    :param eig_val: Given Eigenvalue
    :param points: How much points to calculate for the solution function. This is not required at all.
    :param further: Whether to calculate beyond the boundary parameter. 0 or 1
    :return: parameter list, solution function, 1st derivative of solution function
    """
    Butcher_tableau()
    #   Initializing parameter list with the initial condition.
    x = [init_cond[0]]
    #   Initializing the solution function with the initial condition.
    y1 = [init_cond[1]]
    #   Initializing the 1st derivative of solution function with the initial condition
    y2 = [init_cond[2]]

    #   Main loop to calculate y1 and y2 and increment x with step
    for k in range(int((end_boundary - x[0]) / step)):
        #   Calculating y1
        F(__func_y1, y1[-1], y2[-1], eig_val, x[-1])
        y1.append(y1[-1] + sum([b[i] * F_i[i] for i in range(len(b))]))
        #   Calculating y2
        F(__func_y2, y1[-1], y2[-1], eig_val, x[-1])
        y2.append(y2[-1] + sum([b[i] * F_i[i] for i in range(len(b))]))
        #   Incrementing x with step
        x.append(x[-1] + step)

    #   This is some additional code unused in the program.
    if further:
        index = len(x)
        for k in range(points - index):
            F(__func_y1, y1[-1], y2[-1], eig_val, x[-1])
            y1.append(y1[-1] + sum([b[i] * F_i[i] for i in range(len(b))]))
            F(__func_y2, y1[-1], y2[-1], eig_val, x[-1])
            y2.append(y2[-1] + sum([b[i] * F_i[i] for i in range(len(b))]))
            x.append(x[-1] + step)

    return x, y1, y2


def _get_next_point(eig_val, p_point, p_p_point, x):
    """
        This function calculates the next point of the solution function from the previous
        two points according to the expansion of Sturm-Liouville equation with
        finite difference formula for 1st and 2nd derivative.

    :param eig_val: Given eigenvalue
    :param p_point: 1st previous point of the solution function
    :param p_p_point: 2nd previous point of the solution function
    :param x: Parameter to be evaluated at; the 1st previous x.
    :return: Next point of the solution function
    """
    global p
    global p_d
    global q
    global s

    #   Calculating the coefficients
    beta = (2 * p(x) + step * p_d(x))
    alpha1 = 4 * p(x) - 2 * q(x, eig_val) * step ** 2
    alpha2 = step * p_d(x) - 2 * p(x)
    alpha3 = 2 * s(x) * step ** 2

    return (alpha1 * p_point + alpha2 * p_p_point + alpha3) / beta


def _get_previous_point(eig_val, n_point, n_n_point, x):
    """
        This functions calculates the previous point of the solution function from
        the next two points according to the expansion of the Sturm-Liouville equation
        with the finite differnce formula for 1st and 2nd derivative.

    :param eig_val: Given eigenvalue
    :param n_point: 1st next point of the solution function
    :param n_n_point: 2nd next point of the solution function
    :param x: Parameter to be evaluated at; the 1st next x
    :return: The previous point of the solution function
    """
    global p
    global p_d
    global q
    global s

    #   Calculating the coefficients
    beta = step * p_d(x) - 2 * p(x)
    alpha1 = 4 * p(x) - 2 * q(x, eig_val) * step ** 2
    alpha2 = (2 * p(x) + step * p_d(x))
    alpha3 = 2 * s(x) * step ** 2

    return (-alpha1 * n_point + alpha2 * n_n_point - alpha3) / beta


def _solve_ivp(conditions, eig_val):
    """
        Solving Initial Value Problem with simple Runge-Kutta scheme.

    :param conditions: Conditions of Diff Eq.
    :param eig_val: Given eigenvalue
    :return: Parameter and the solution function
    """
    X, Y, Y_d = _runge_kutta(conditions, conditions[0] + step, eig_val, step, 0)
    return X, Y


def _solve_bvp(conditions, eig_val):
    """
        Solving Boundary Value Problem with Shooting method.

    :param conditions: Conditions of Diff Eq.
    :param eig_val: Given eigenvalue
    :return: Parameter and the solution function
    """

    #   Settnig the tolerance as zero
    zero = 1e-5
    #   Initial guess for the secant method to run at
    guess = 0.1

    # Isolating the initial condition from conditions
    init_cond = conditions[0:2]
    # Isolating the boundary condition from conditions
    bound_cond = conditions[2:4]

    #   Assuming the initial value of the 1st derivative of the solution function
    #   as guess.
    init_cond.append(guess)

    #   Approximating a second guess
    guess2 = guess + step
    guess3 = 0

    #   The loop for executing secant method
    while abs(guess2 - guess) >= zero:
        init_cond[2] = guess
        x, y1, y2 = _runge_kutta(init_cond, bound_cond[0], eig_val, step, 0)
        y1_boundary_f = y1[-1] - bound_cond[1]
        init_cond[2] = guess2
        x, y1, y2 = _runge_kutta(init_cond, bound_cond[0], eig_val, step, 0)
        y1_boundary2_f = y1[-1] - bound_cond[1]

        dnmntr = (y1_boundary2_f - y1_boundary_f)
        guess3 = guess2 - y1_boundary2_f * (guess2 - guess) / dnmntr

        guess = guess2
        guess2 = guess3

    #   Found the initial value of the 1st derivative of the solution function
    #   to be guess3
    init_cond[2] = guess3

    #   Executing Runge-Kutta for proper initial conditions
    X, Y, Y_d = _runge_kutta(init_cond, bound_cond[0], eig_val, step, 0)

    return X, Y


def get_func(eig_val, plot_range):
    """
        The main function of the program to get the solution function for
        given conditions and eigenvalue.

        First it checks whether the problem is an IVP or a BVP.

        Then it calls _solve_ivp or _solve_bvp accordingly and takes only two values
        of the solution function.

        With those two values it calculates other values of solution function for the
        parameter between given plot range with _get_next_point and _get_previous_point

    :param eig_val: Given eigenvalue
    :param plot_range: Parameter range for the solution function to be evaluated at
    :return: Parameter, Solution function
    """
    global conditions

    #   Initializing X and Y (The solution function)
    X = 0
    Y = 0

    #   Checking for the problem type and solve accordingly
    if PROB == 'IVP':
        X, Y = _solve_ivp(conditions, eig_val)

    elif PROB == 'BVP':
        X, Y = _solve_bvp(conditions, eig_val)

    #   Taking the first two values of X acquired from solving particular
    #   differential equation problem
    X = X[0:2]
    Y = Y[0:2]

    #   X_p and Y_p are declared for holding the values of the solution function
    #   where parameter x is decremented to go through the left of the plot to
    #   get the values between plot range.
    #   The more the index get higher, the less the x becomes. Thus, X_p and Y_p
    #   are taken as the reverse of X and Y
    X_p = copy.deepcopy([X[1], X[0]])
    Y_p = copy.deepcopy([Y[1], Y[0]])
    #   Obtaining Y_p (and X_p) with _get_previous_point
    while X_p[-1] > plot_range[0]:
        Y_p.append(_get_previous_point(eig_val, Y_p[-1], Y_p[-2], X_p[-1]))
        X_p.append(X_p[-1] - step)

    #   Reversing X_p and Y_p to be able to merge them with X and Y later.
    #   Think of it as sorting X_p from low value to high value
    X_p.reverse()
    Y_p.reverse()

    #   Obtaining Y (and X) with _get_next_point
    while X[-1] < plot_range[1]:
        Y.append(_get_next_point(eig_val, Y[-1], Y[-2], X[-1]))
        X.append(X[-1] + step)

    #   Merging X with X_p and Y with Y_p in a sorted order and returning
    return X_p[2:] + X, Y_p[2:] + Y


def legendre_equation(eig_val):
    """
        This function defines the functions required for the Sturm-Liouville equation
        to be Legendre equation and asks user for the conditions in order to solve it.


    :param eig_val: Given eigenvalue
    :return: Proper eigenvalues for the general Sturm-Liouville equation
    """
    global p
    global p_d
    global q
    global s
    global conditions
    global PROB

    p = lambda x: (1 - x ** 2)
    p_d = lambda x: -2 * x
    r = lambda x: 0
    w = lambda x: 1
    q = lambda x, eig_value: -r(x) + eig_value * w(x)
    s = lambda x: 0
    eigenvalue = lambda x: x * (x + 1)
    eigenvalues = [eigenvalue(i) for i in eig_val]
    __set_conditions()

    return eigenvalues


def bessel_equation(eig_val):
    """
        This function defines the functions required for the Sturm-Liouville equation
        to be Bessel equation and asks user for the conditions in order to solve it.


    :param eig_val: Given eigenvalue
    :return: Proper eigenvalues for the general Sturm-Liouville equation
    """
    global p
    global p_d
    global q
    global s
    global conditions
    global PROB

    p = lambda x: x
    p_d = lambda x: 1
    r = lambda x: -x
    w = lambda x: -1 / x
    q = lambda x, eig_value: -r(x) + eig_value * w(x)
    s = lambda x: 0
    eigenvalue = lambda x: x ** 2
    eigenvalues = [eigenvalue(i) for i in eig_val]
    __set_conditions()

    return eigenvalues


def custom_equation(eig_val):
    """
        This function asks user for the required functions and the conditions
        of a custom Sturm-Liouville equation.


    :param eig_val: Given eigenvalue
    :return: Proper eigenvalues for the general Sturm-Liouville equation
    """
    global p
    global p_d
    global q
    global s
    global conditions

    x, y, z = sp.symbols("x y z")
    p = sp.sympify(input("Enter p(x) : "))
    p = sp.lambdify(x, p, 'numpy')
    p_d = sp.sympify(input("Enter p'(x) : "))
    p_d = sp.lambdify(x, p_d, 'numpy')
    r = sp.sympify(input("Enter r(x) : "))
    r = sp.lambdify(x, r, 'numpy')
    w = sp.sympify(input("Enter w(x) : "))
    w = sp.lambdify(x, w, 'numpy')
    s = sp.sympify(input("Enter s(x) : "))
    s = sp.lambdify(x, s, 'numpy')
    q = lambda x, eig_value: -r(x) + eig_value * w(x)
    eigenvalue_ev = sp.sympify(input("Enter eigenvalue evaluator : "))
    eigenvalue_ev = sp.lambdify(x, eigenvalue_ev, 'numpy')
    eigenvalues = [eigenvalue_ev(i) for i in eig_val]
    __set_conditions()

    return eigenvalues


if __name__ == '__main__':
    eq = int(input("Which function do you want to plot \n "
                   "\t[1] Legendre Equation\n"
                   "\t[2] Bessel Equation\n"
                   "\t[3] Custom Equation\n\n"
                   "\tSelect (1/2/3): "))
    print('____________________________________________________')
    #   Asking user for the eigenvalues
    eig_val = list(map(float, [i for i in input("Enter some eigenvalues "
                                                "(with comma to separate): ").split(',')]))
    #   Asking user for the plot range
    plot_range = list(map(float, [i for i in input('Enter plot range '
                                                   '(two values of x separated with comma): ').split(',')]))
    #   Sort the plot range
    plot_range.sort()

    eigenvalues = 0
    if eq == 1:
        eigenvalues = legendre_equation(eig_val)
    elif eq == 2:
        eigenvalues = bessel_equation(eig_val)
    elif eq == 3:
        eigenvalues = custom_equation(eig_val)

    #   Initializing plot
    fig, ax = plt.subplots()
    plt.xlabel("x")
    plt.ylabel("Y")
    plt.xlim(plot_range[0], plot_range[1])
    # plt.ylim(-1.5, 1.5)
    #   Solving and plotting for the given eigenvalues
    i = 0
    for eigenvalue in eigenvalues:
        X, Y = get_func(eigenvalue, plot_range)
        ax.plot(X, Y, label='Y {}'.format(eig_val[i]), marker='o', markersize=1)
        i += 1
        fig.legend()
        fig.show()
    plt.show()
