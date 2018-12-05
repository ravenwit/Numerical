import sys
import numpy as np
import copy as cp
from matplotlib import pyplot as plt
import matplotlib.animation as anim
import mpl_toolkits.mplot3d.axes3d as p3
import RungeKutta as rk

'''
trajectory = np.array([[[0, 0, 0]], [[2,2,2]]])
np.concatenate((trajectory[1], np.array([[1,2,3]])), axis=0)
trajectory = np.array([trajectory[0], np.concatenate((trajectory[1], np.array([[65,76,89]])), axis=0)])
array([[2, 2, 2],
       [1, 2, 3]])
'''

h = 0.00001
points = 1000

c = 137
q = [1, 1]
m = [1836, 1]

time = []
trajectory = np.array([[], []])
velocity = np.array([])
acceleration = np.array([[], []])
radiation = np.array([[], []])
hamil = np.array([[], []])


def _init_trajectory():
    global trajectory
    global velocity
    global acceleration
    global radiation
    global time

    shaped_zeroes = np.zeros(points * 3).reshape(points, 3)
    trajectory = np.array([shaped_zeroes, shaped_zeroes])
    velocity = np.array([shaped_zeroes, shaped_zeroes])
    acceleration = np.array([shaped_zeroes, shaped_zeroes])
    radiation = np.array([np.zeros(points), np.zeros(points)])

    # radius = 0.005
    # radius = 0.004
    # radius = 0.005
    radius = 0.05
    # theta = np.pi / 4
    theta = np.pi / 4
    phi = np.pi / 4
    # phi = np.pi / 2
    x_radius = radius * np.sin(theta) * np.cos(phi)
    y_radius = radius * np.sin(theta) * np.sin(phi)
    z_radius = radius * np.cos(phi)

    trajectory[0][0] = [0, 0, 0]
    trajectory[1][0] = [x_radius, y_radius, z_radius]

    velocity[0][0] = [0, 0, 0]
    # velocity[1][0] = [0.001 * c, -0.008 * c, -0.001 * c]
    # velocity[1][0] = [-0.08 * c, -0.0008 * c, 0.00001 * c]
    # velocity[1][0] = [-0.03 * c, -0.0003 * c, 0.00001 * c]
    # velocity[1][0] = [0.04 * c, -0.0008 * c, 0.00001 * c]
    # velocity[1][0] = [0.1 * c, 0.1 * c, 0.1 * c]
    velocity[1][0] = [0.01 * c, -0.01 * c, -0.01 * c]
    # velocity[1][0] = [0.001 * c, -0.008 * c, -0.001 * c]

    acceleration[0][0] = [0, 0, 0]
    acceleration[1][0] = [0, 0, 0]

    time = [0]


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


def _retarded_index(body, current_index, position, t):
    global h
    # zero = 1e-5

    guess_index = current_index - 1
    guess = time[guess_index]
    guess2_index = current_index - 2
    guess2 = time[guess2_index]

    while abs(guess2_index - guess_index) != 0:
        f1 = c * (t - guess) - np.linalg.norm(position - trajectory[int(not body)][guess_index])
        f2 = c * (t - guess2) - np.linalg.norm(position - trajectory[int(not body)][guess2_index])

        dnmntr = (f2 - f1)
        guess3 = guess2 - f2 * (guess2 - guess) / dnmntr

        guess = guess2
        guess_index = guess2_index
        guess2 = guess3
        guess2_index = int(guess3 * (1 / h))
        if guess2_index <= 0:
            guess2_index = 0
            break

    return guess2_index


def EBfield(body, index, position, t):
    global h
    global trajectory
    global velocity
    global acceleration
    global time
    global points
    E = 0
    B = 0
    if index < 2:
        R = position - trajectory[int(not body)][index]
        r = np.linalg.norm(R)
        R_unit = R / r
        E = q[int(not body)] / r ** 2
        E = E * R_unit
        B = np.cross(R_unit, E) / c

    else:
        ret_index = _retarded_index(body, index, position, t)

        R = position - trajectory[int(not body)][ret_index]
        r = np.linalg.norm(R)
        R_unit = R / r

        V = velocity[int(not body)][ret_index]
        v = np.linalg.norm(V)
        beta = V / c
        gamma = 1 / np.sqrt(1 - v ** 2 / c ** 2)
        # V_h = velocity[int(not body)][ret_index - 1]
        # V_dot = (V - V_h) / h
        V_dot = acceleration[int(not body)][ret_index]
        beta_dot = V_dot / c

        VR = np.dot(V, R)

        # factor1 = R_unit - beta
        # factor2 = (1 - np.dot(R_unit, beta)) ** 3
        # factor3 = (gamma ** 2) * factor2 * r ** 2
        # factor4 = factor1 / factor3
        # factor5 = np.cross(R_unit, np.cross(factor1, beta_dot))
        # factor6 = c * factor2 * r
        # factor7 = factor5 / factor6
        # factor8 = factor4 + factor7
        # # factor8 = factor4
        #
        # E = q[int(not body)] * factor8

        factor1 = 1 / (r - VR / c) ** 3
        factor2 = R - r * V / c
        factor3 = 1 - (v / c) ** 2
        factor4 = factor3 * factor2
        factor5 = V_dot / c ** 2
        factor6 = np.cross(factor2, factor5)
        factor7 = np.cross(R, factor6)
        factor8 = factor4 + factor7
        # factor8 = factor4
        factor9 = factor1 * factor8

        E = q[int(not body)] * factor9

        #
        B = np.cross(R_unit, E)

    return E, B


def LorentzForce(body, index, position, velocity, t):
    E, B = EBfield(body, index, position, t)
    V = velocity

    # if index>4: return q[body] * (E + np.cross(V, B)) - (((acceleration[body][index]-acceleration[body][index-1])/h)*q[body]**2)/(6*np.pi*c)
    # else: return q[body] * (E + np.cross(V, B))
    # return np.array([3260,5230,-630])
    return q[body] * (E + np.cross(V, B))


def __velocity(position, velocity, t, **kwargs):
    return velocity


def __acceleration(position, velocity, t, **kwargs):
    body = 0
    index = 0
    if kwargs is not None:
        for key, value in kwargs.items():
            if key == "body":
                body = value
            elif key == "index":
                index = value
    else:
        sys.exit()

    beta = np.linalg.norm(velocity) / c
    # factor1 = np.sqrt(1 - beta ** 2)
    # factor2 = 1 + beta ** 2 / (1 - beta ** 2)
    # factor3 = factor2 * m[body]
    # factor4 = factor1 / factor3

    gamma = 1 / np.sqrt(1 - beta ** 2)
    F = LorentzForce(body, index, position, velocity, t)
    factor1 = 1 / (m[body] * gamma)
    factor2 = np.dot(velocity, F)
    factor3 = factor2 * velocity
    factor4 = factor3 / c ** 2
    factor5 = F - factor4
    acc = factor1 * factor5

    return acc


def _ode_solve(body, index):
    global h
    global trajectory
    global velocity
    global time

    last_position = trajectory[body][index]
    last_velocity = velocity[body][index]

    last_time = time[index]

    n_t, n_position, n_velocity = rk.get_single_point(last_position, last_velocity, last_time, __velocity,
                                                      __acceleration, h, 1, body=body, index=index)
    return n_position, n_velocity



# ax.set_xlim3d([0.0, 0.01])


def get_trajectory_iter(index):
    global h
    global trajectory
    global velocity
    global acceleration
    global time

    if index == 0:
        _init_trajectory()
    else:

        # For proton
        n_position0, n_velocity0 = _ode_solve(0, index - 1)

        # print(n_position0)

        trajectory[0][index] = n_position0
        velocity[0][index] = n_velocity0

        # For electron
        n_position1, n_velocity1 = _ode_solve(1, index - 1)

        print(n_position1)

        trajectory[1][index] = n_position1
        velocity[1][index] = n_velocity1

        time.append(time[-1] + h)

        acceleration[0][index] = __acceleration(trajectory[0][index], velocity[0][index], time[index],
                                                body=0, index=index)
        acceleration[1][index] = __acceleration(trajectory[1][index], velocity[1][index], time[index],
                                                body=1, index=index)


def get_trajectory_anim(index):
    global h
    global trajectory
    global velocity
    global acceleration
    global time

    if index == 0:
        _init_trajectory()
    else:
        # For proton
        n_position0, n_velocity0 = _ode_solve(0, index - 1)

        print(n_position0)

        trajectory[0][index] = n_position0
        velocity[0][index] = n_velocity0

        # For electron
        n_position1, n_velocity1 = _ode_solve(1, index - 1)

        print(n_position1)

        trajectory[1][index] = n_position1
        velocity[1][index] = n_velocity1

        time.append(time[-1] + h)
        index += 1
        # ax.clear()
        # ax.plot(trajectory[0][:index, 0], trajectory[0][:index, 1], trajectory[0][:index, 2], c='blue', marker='o',
        #         markersize=9)
        # ax.plot(trajectory[1][:index, 0], trajectory[1][:index, 1], trajectory[1][:index, 2], c='red', marker='o',
        #         linestyle='none',
        #         markersize=3)
        # ax.plot(np.array([trajectory[1][:index, 0][-1]]), np.array([trajectory[1][:index, 1][-1]]),
        #         np.array([trajectory[1][:index, 2][-1]]), c='green', marker='>',
        #         markersize=9)
        # if index: plt.savefig('twobody/trajec1/figure{}'.format(index))


def get_trajectory():
    for i in range(points):
        get_trajectory_anim(i)


def plotEnergy():
    global time

    timeE = cp.deepcopy(time)
    E = []
    for i in range(len(timeE)):
        V = velocity[0][i]
        v = np.linalg.norm(V)
        beta = v / c
        gamma = 1 / np.sqrt(1 - beta ** 2)
        m0 = m[0] * gamma
        E0 = m0 * c ** 2

        V = velocity[1][i]
        v = np.linalg.norm(V)
        beta = v / c
        gamma = 1 / np.sqrt(1 - beta ** 2)
        m1 = m[1] * gamma
        E1 = m1 * c ** 2
        E.append((E0 + E1))
        # if i>0:
        #     E.append(E[0]-(E0 + E1))
        # else:
        #     E.append((E0 + E1))

    timeE = timeE[3:]
    maxE = max(E)
    E = np.array(E) / maxE
    E = np.array(E[3:])
    maxE = max(E)
    minE = min(E)
    E = np.array(list(map(lambda x: translate(x, minE, maxE, 0, 1000), E)))

    fig1, ax1 = plt.subplots()
    plt.title("Time Energy")
    plt.xlabel("Time")
    # plt.ylim(99.99, 100)
    plt.ylabel("Energy")
    ax1.plot(timeE, E)
    fig1.legend()
    fig1.show()
    plt.savefig("twobody/Energy.png")


def plotTrajectory():
    fig = plt.figure()
    # ax = p3.Axes3D(fig)
    ax = fig.gca(projection='3d')

    # ax.set_title("Trajectories", fontsize=18)
    ax.set_xlabel("X Axis", fontsize=16)
    ax.set_ylabel("Y Axis", fontsize=16)
    ax.set_zlabel("Z Axis", fontsize=16)
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    ax.set_xticklabels([])

    ax.plot(trajectory[0][:, 0], trajectory[0][:, 1], trajectory[0][:, 2], c='blue', marker='o',
            markersize=9)
    ax.plot(trajectory[1][:, 0], trajectory[1][:, 1], trajectory[1][:, 2], c='red', marker='o',
            markersize=3)
    ax.plot(np.array([trajectory[0][:, 0][-1]]), np.array([trajectory[0][:, 1][-1]]),
            np.array([trajectory[0][:, 2][-1]]), c='green', marker='>',
            markersize=9)
    ax.plot(np.array([trajectory[1][:, 0][-1]]), np.array([trajectory[1][:, 1][-1]]),
            np.array([trajectory[1][:, 2][-1]]), c='magenta', marker='>',
            markersize=9)

    fig.show()
    plt.savefig("twobody/trajectory.png")


def plotTrajectoryR():
    # ax.plot(trajectory[0][:, 0], trajectory[0][:, 1], trajectory[0][:, 2], c='blue', marker='o',
    #         markersize=9)
    # ax.plot(trajectory[1][:, 0], trajectory[1][:, 1], trajectory[1][:, 2], c='red', marker='o',
    #         markersize=3)
    ax.plot(np.array([trajectory[0][:, 0][-1]]), np.array([trajectory[0][:, 1][-1]]),
            np.array([trajectory[0][:, 2][-1]]), c='blue', marker='o',
            markersize=15)
    ax.plot(np.array([trajectory[1][:, 0][-1]]), np.array([trajectory[1][:, 1][-1]]),
            np.array([trajectory[1][:, 2][-1]]), c='red', marker='o',
            markersize=15)

    fig.show()
    plt.savefig("twobody/trajectory.png")


def plotVelocity(body):
    timeE = cp.deepcopy(time)
    timeE = timeE[2:]
    velocityE = cp.deepcopy(velocity)
    velocityE = velocityE[:, 2:]
    # velocityE = velocityE

    fig1, ax1 = plt.subplots()
    # plt.title("Charge: {} ; Mass: {}".format(q[body], m[body]))
    plt.xlabel("Time")
    plt.ylabel("Velocity")
    ax1.plot(timeE, [np.linalg.norm(velocityE[body][i]) for i in range(len(timeE))])
    fig1.show()
    plt.savefig("twobody/Velocity.png")


def plotVelocityRadius(body):
    RV2 = []
    for i in range(len(time)):
        rad = np.linalg.norm(trajectory[body][i])
        v2 = np.linalg.norm(trajectory[body][i]) ** 2
        RV2.append(rad * v2)

    fig1, ax1 = plt.subplots()
    plt.title("Charge: {} ; Mass: {}".format(q[body], m[body]))
    plt.xlabel("Time")
    plt.ylabel("Radius * Velocity^2")
    ax1.plot(time, RV2)
    fig1.show()
    plt.savefig("twobody/VeloctyRadius.png")


def plotRadius(body):
    RV2 = []
    for i in range(len(time)):
        rad = np.linalg.norm(trajectory[int(not body)][i] - trajectory[body][i])
        RV2.append(rad)

    fig1, ax1 = plt.subplots()
    # plt.title("Charge: {} ; Mass: {}".format(q[body], m[body]))
    plt.xlabel("Time")
    plt.ylabel("Radius")
    ax1.plot(time, RV2)
    fig1.show()
    plt.savefig("twobody/radius.png")


def plotMomentum():
    global time

    M = []

    for i in range(len(time)):
        V = velocity[0][i]
        v = np.linalg.norm(V)
        beta = v / c
        gamma = 1 / np.sqrt(1 - beta ** 2)
        m0 = m[0] * gamma
        M0 = m0 * v

        V = velocity[1][i]
        v = np.linalg.norm(V)
        beta = v / c
        gamma = 1 / np.sqrt(1 - beta ** 2)
        m1 = m[1] * gamma
        M1 = m1 * v
        M.append((M0 + M1))
        # if i>0:
        #     M.append(M[0]-(M0 + M1))
        # else:
        #     M.append((M0 + M1))

    maxM = max(M)
    M = np.array(M) / maxM
    maxM = max(M)
    minM = min(M)
    M = np.array(list(map(lambda x: translate(x, minM, maxM, 0, 1000), M)))

    fig1, ax1 = plt.subplots()
    plt.title("Time Momentum")
    plt.xlabel("Time")
    # plt.ylim(99.99, 100)
    plt.ylabel("Momentum")
    ax1.plot(time, M)
    fig1.legend()
    fig1.show()
    plt.savefig("twobody/momentum.png")


def plotRadiation(body):
    M = []

    for i in range(len(time)):
        V = velocity[0][i]
        v = np.linalg.norm(V)
        beta = v / c
        gamma = 1 / np.sqrt(1 - beta ** 2)
        m0 = m[0] * gamma
        M0 = m0 * v

        V = velocity[1][i]
        v = np.linalg.norm(V)
        beta = v / c
        gamma = 1 / np.sqrt(1 - beta ** 2)
        m1 = m[1] * gamma
        M1 = m1 * v
        M.append((M0 + M1))
        radiation[body][i] = (2 / 3) * q[body] ** 2 / (m[body] ** 2 * c ** 3) * ((M[i] - M[i - 1] / h) ** 2)
        # if i>0:
        #     M.append(M[0]-(M0 + M1))
        # else:
        #     M.append((M0 + M1))

    fig1, ax1 = plt.subplots()
    # plt.title("Charge: {} ; Mass: {}".format(q[body], m[body]))
    plt.xlabel("Time (in sec.)")
    plt.ylabel("Radiation Power (in Ha/sec)")
    ax1.plot(time, radiation[body])
    fig1.show()
    plt.savefig("twobody/radiation.png")


def plotYZ(body):
    fig1, ax1 = plt.subplots()
    plt.title("Charge: {} ; Mass: {}".format(q[body], m[body]))
    plt.xlabel("Z")
    plt.ylabel("X")
    ax1.plot(time, trajectory[body][:, 1])
    fig1.show()


def plotVE(body):
    global time

    E = []
    for i in range(len(time)):
        V = velocity[body][i]
        v = np.linalg.norm(V)
        beta = v / c
        gamma = 1 / np.sqrt(1 - beta ** 2)
        m0 = m[body] * gamma
        E0 = m0 * c ** 2

        E.append(E0)
        # if i>0:
        #     E.append(E[0]-(E0 + E1))
        # else:
        #     E.append((E0 + E1))

    maxE = max(E)
    minE = min(E)
    E = np.array(list(map(lambda x: translate(x, minE, maxE, 0, 1000), E)))

    fig1, ax1 = plt.subplots()
    plt.title("Charge: {} ; Mass: {}".format(q[body], m[body]))
    plt.xlabel("Velocity")
    plt.ylabel("Energy")
    ax1.plot([np.linalg.norm(V) for V in velocity[body]], E)
    fig1.show()
    plt.savefig("twobody/velocityEnergy.png")


if __name__ == '__main__':
    # ani = anim.FuncAnimation(fig, get_trajectory_anim, points, interval=1, blit=False)
    get_trajectory()
    # plotVE(1)
    # plotEnergy()
    # plotMomentum()
    plotTrajectory()
    # plotTrajectoryR()
    # plotRadiation(0)
    # plotRadiation(1)
    # plotVelocity(0)
    plotVelocity(1)
    plotRadius(1)
    # plotVelocityRadius(1)
    # plotYZ(1)
    print('Upadted')
    plt.show()
