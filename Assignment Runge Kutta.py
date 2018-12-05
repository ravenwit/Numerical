"""
    This code comes with an article on 'obtaining the formulae for
    12th order Runge-Kutta method'
    titled
    "Developing an algorithm for 12th Order Runge-Kutta Method"

    The theoretical background to understand the inner mechanism and
    numerical implementation of the following script can be found in
    the article.
"""



##
#   Importing necessary libraries
##
import numpy as np
import matplotlib.pyplot as plt
#####################################################################


##
#   Declaring two enum to indicate which order to use for Runge-Kutta
#   4th order is indicated as RK4 and 12th order is indicated as RK12
##
RK4 = 1
RK12 = 2
#####################################################################


##
#   Initializing some universal list to put values into them while program
#   runs.
#
#   F_i holds the value for F_i
#   a holds the value of A matrix elements, a_ij
#   b holds the value of weight coefficients, b_i
#   c holds the value of time increment within a single step, c_i
##
F_i = []
a = []
b = []
c = []
#####################################################################


##
#  Butcher Tableau
##
def Butcher_tableau(method):
    """
        This function just inserts the value of Butcher Matrix to its
        elements a, b, c.

        It takes an argument to identify which order (4 or 12) is being
        used and puts the according value to a, b, c which are pre-written.

    :param method: RK4 or RK12
    :return: Returns the arrays a, b, c. Which isnt necessary. a, b and c are global
             variable.
    """

    # Referring to the global versions of a, b, c and F_i
    global a
    global b
    global c
    global F_i

    # If method equals 1 then the order is 4
    if method == 1:
        c = np.array([0, 1 / 2, 1 / 2, 1])
        a = np.array([
            [0, 0, 0, 0],
            [1 / 2, 0, 0, 0],
            [0, 1 / 2, 0, 0],
            [0, 0, 1, 0]
        ])
        b = np.array([1 / 6, 1 / 3, 1 / 3, 1 / 6])

    #  If method is 2 then the order is 12
    elif method == 2:
        c = np.array([0,0.2,5/9,5/6,1/3,1,67801993/100920496,44758598/155021585,9/16,5/6,172135849/181636255,36851109/672327007,20592542/242584693,
                      30120495 / 113415896,1/2,210568577/286712394,96134905/105052617,172135849/181636255,5/6,44758598/155021585,67801993/100920496,1/3,5/9,1/5,1])
        a=np.array([
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            [1/5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
            [-35/162, 125/162, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [5/24, 0, 5/8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [29/150, 0, 11/50, -2/25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1/10, 0, 0, 2/5, 1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [37074941/358681667, 0, 0, 33585619/270735842, 257720006/533392767, -44020349/1135920344, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [25230124/203405979, 0, 0, 0, 98101881/451976942, 65957727/4798468366, -6070937/91831493, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [24737624/270423075, 0, 0, 0, 0, -3736153/686353106, 19485332/286247261, 50525157/123716602, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [60247834/676931571, 0, 0, 0, 0, 4275247/855856941, 36672659/92161292, 16731755/39099261, -38981915/450596697, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [27288061/392584475, 0, 0, 0, 0, 365006343/2826287155, 141373675/92356644, 163264999/282526613, -51629527/54272901, 142621687/349326099, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [32856654/738581809, 0, 0, 0, 0, -14162705/3722356397, 25582922/2391931811, 51119620/2438724161, -25704425/1102503257, 4688654/1780957031, 3117485/988194642, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [10458484/537465835, 0, 0, 0, 0, 0, 0, 0, 296227/4365826774, -383576/8924609019, 146141/8286564037, 12646113/193405084, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [31431986/151965127, 0, 0, 0, 0, 0, 0, 0, 8107993/486102169, -20936474/2380493097, 11794237/3402097500, -152083614/176581783, 102389113/112682442, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [11253958/551864565, 0, 0, 0, 0, 0, 0, 0, 42975146/494268647, -4474701/233483414, 6268097/956043048, 15766372/159663323, 117183268/21888493765, 28144139/93450007, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [89861716/393422115, 0, 0, 0, 0, 0, 0, 0, -190920249/382830190, 37407972/277422485, -15879177/409829375, -24748163/19414396, 49033497/34070828,-135425641/632808015, 100437343/104818503, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [329422004/164527983, 0, 0, 0, 0, 0, 0, 0, 159202720/77020477, 152179118/243885337, -35567233/769381099, -466749240/52741619, 161194222/20819195, -453082123/770078291,-26419211/23869100,-79543028/85573473, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [121929017/38856942, 0, 0, 0, 0, 365006343/2826287155, 141373675/92356644, 163264999/282526613, 439144222/81009727, 13186261/56948547, 25264555/332737891, -401683667/32464540, 2492304541/252908789, 47761777/555943912, -75613337/13377144, -248933559/128117530, -36061621/280957460, 0, 0, 0, 0, 0, 0, 0, 0],
            [155993079/112744303, 0, 0, 0, 0, 4275247/855856941, 36672659/92161292,16731755/39099261, -164168531/125993596, 855875317/1294246651, -43179215/298694538, -509694778/73171449, 235461367/35364726, -82607643/49466432,119839851/58058089, -108798317/161243854,-1216233/1051933279, -11104923/2041128862, 0, 0, 0, 0, 0, 0, 0],
            [30575399/32142801, 0, 0, 0, 98101881/451976942, 65957727/4798468366, -6070937/91831493, 0, 58237611/382433426, -87760066/259844263, -9942595/515625276, -123996494/33670977, 101639181/32144170, -92687960/250195241, -4868661/94541843, -870114/1048803305,240782/86054719017,3179847/75963145,24618431/88211465, 0, 0, 0, 0, 0, 0],
            [37074941 / 358681667, 0, 0, 33585619 / 270735842, 257720006 / 533392767, -44020349 / 1135920344, 0,-120527899 / 274980832, 0, -58596715 / 268009592, -47102449 / 1508075769, 0, 0, 0, 0, 0, 0,47102449 / 1508075769, 58596715 / 268009592, 120527899 / 274980832, 0, 0, 0, 0, 0],
            [29/150, 0, 11/50, -2/25, 0,0, 14340957/145703507, -31420077/159971156, 0, 103497843/237131315, 22131167/339115870, 0, 0, 0, 0, 0, 0, -22131167/339115870, -103497843/237131315,31420077/159971156, -14340957/145703507, 0, 0, 0, 0],
            [-35/162, 125/162, 0, 0, -2/3, 0, -82382086/210859561, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 82382086/210859561,2/3, 0, 0, 0],
            [1/5, 0, -40/243, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 40/243, 0, 0],
            [157084639/106730534, 63/80, 91/216, 0,7/24, 0, 100146079/287280186, 70568068/307486745,129515617/22367050,72807936/173937191, 51739211/168509742, -123960914/26447765 ,207695632/66235459, 222633283/158870770,-234969570/42495271, -108932901/127684936,28937867/279388356, -46373597/330121299, -72807936/173937191, -70568068/307486745, -100146079/287280186, -7/24, -91/216, -63/80, 0]
        ])
        b=np.array([1/42,3/128,1/32,0,1/24,0,1/20,1/20,0,1,1/14,0,38532146/278385263,
                    119008733/551291285,128/525,119008733/551291285,38532146/278385263,-1/14,
                    -1/10,-1/20,1/20,-1/24,-1/32,-3/128,1/42])

    # Initializing F_i with zeroes. The length of F_i is equal to the length of b.
    # Thus initializing initializing with len(b) number of zeroes
    F_i = np.array(np.zeros(len(b)))
    return a, b, c


def __func_y1(y1, y2, t):
    """
        The y1, which is actually theta in the theory.
        The derivative of y1 returns a function y2 which is just the derivative of y1
        namely angular velocity.

        The need for y1 and t in this function is that we will later use a general
        indicator for both __func_y1 and __finc_y2. Therefore, their arguments will
        have to same.

    :param y1: y1, Theta at time t
    :param y2: y2, Angular velocity at time t
    :param t: Time
    :return: y2, Angular velocity at time t
    """
    return y2


def __func_y2(y1, y2, t):
    """
        This fucntion is the right hand side of the equation of 'Force Driven Chaotic
        Pendulum' if the left side is only the 2nd order differential, derivative of
        y2.

        q, w and omega0 are parameters for the physical pendulum. They are fined tuned
        here.

    :param y1: y1, Theta at time t
    :param y2: y2, Angular velocity at time t
    :param t: Time
    :return: f(t, y1, y2)
    """
    q = 1.15
    w = 0.9
    omega0 = 2.0 / 3.0
    return -q * y2 - np.sin(y1) + w * np.cos(omega0 * t)


def F(func, y1, y2, t, step):
    """
        Calculating the value of F_i for a single time step.

        We are actually solving two differential equations here to obtain a solution
        for the 2nd order differential.

    :param func: Which function to use, __func_y1 or __func_y2.
    :param y1: y1, Theta at time t
    :param y2: y2, Angular velocity at time t
    :param t: Time
    :param step: h, one time step
    :return: returns nothing. But update the global variable F_i for latter use.
    """
    # Referring to the global version of F_i
    global F_i
    # Calculating according to the theory. The article has a detailed explanation.
    for i in range(len(b)):
        F_i[i] = step * func((y1 + sum([a[i][j] * F_i[j] for j in range(len(a))])),
                               (y2 + sum([a[i][j] * F_i[j] for j in range(len(a))])),
                               t + step * c[i])


def get_func(y10, y20, t0, step, points, method):
    """
        This function gets the function y, the solution to the differential equation.

    :param y10: The value for y1 at t0 (Initial value)
    :param y20: The value for y2 at t0 (Initial value)
    :param t0: The initial value for time t
    :param step: h, time step
    :param points: integer; How many points to calculate for solution function y
    :param method: RK$ or RK12; Which order to use for Runge-Kutta
    :return:
    """
    #   Initializing the values a, b, c with the proper values for order indicated by the
    #   parameter method
    Butcher_tableau(method)
    #   Putting the initial value for t at the first of the list
    t = [t0]
    #   Putting the initial value for y1 at the first of the list
    y1 = [y10]
    #   Putting the initial value for y2 at the first of the list
    y2 = [y20]
    #   Calculating the solution function with time increment h=step
    for k in range(points):
        #   Initializes F_i with proper values for y1
        F(__func_y1, y1[-1], y2[-1], t[-1], step)
        #   Calculating the next y1 from previous y1
        y1.append(y1[-1] + sum([b[i] * F_i[i] for i in range(len(b))]))
        #   Initializes F_i with proper values for y2
        F(__func_y2, y1[-1], y2[-1], t[-1], step)
        #   Calculating the next y2 from previous y1
        y2.append(y2[-1] + sum([b[i] * F_i[i] for i in range(len(b))]))
        #   Incrementing time
        t.append(t[-1] + step)

    return t, y1, y2


if __name__ == '__main__':
    #   Getting a solution space with y10=0, y20=1, t0=0, step=0.1, points=100, method=RK4)
    t, y1, y2 = get_func(0, 1, 0, 0.1, 100, RK4)
    #   Getting a solution space with y10=0, y20=1, t0=0, step=0.1, points=100, method=RK12)
    t_m2, y1_m2, y2_m2 = get_func(0, 1, 0, 0.1, 100, RK12)
    # plt.plot(t, y1, color='#591919', linewidth=2, label='Theta(4)', marker='o', markersize=1)
    # plt.plot(t, y2, color='#38aac4', linewidth=2, label='Angular velocity(4)', marker='+', markersize=1)
    # plt.plot(t_m2, y1_m2, color='#34bc54', linewidth=2, label='Theta(12)', marker='o', markersize=1)
    # plt.plot(t_m2, y2_m2, color='#b8d13c', linewidth=2, label='Angular velocity(12)', marker='+', markersize=1)
    #   Plotting the graph (y1 vs. y2)
    plt.xlabel('Theta(y1)')
    plt.ylabel('Angular velocity(y2)')
    fig = plt.plot(y1, y2, label='RK4', marker='o', markersize=5)
    fig2 = plt.plot(y1_m2, y2_m2, label='RK12', marker='o', markersize=5)
    plt.legend()
    plt.show()

