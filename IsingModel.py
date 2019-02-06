#!/usr/bin/env python

"""
    Importing necessary libraries
"""

from __future__ import division
import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import cm
import matplotlib.animation as anim


class Ising():
    """
        Definitions
    """

    #   Defining Plot
    f = plt.figure(figsize=(10, 10), dpi=80)
    ax = f.gca()

    #   Physical Constants
    N_t = 88  # Number  of temperature points
    N = 16  # Size of lattice (N x N)
    T = np.linspace(1.53, 3.28, N_t)  # Temperature values
    E = np.zeros(N_t)  # Energy values
    M = np.zeros(N_t)  # Magnetization values
    C = np.zeros(N_t)  # Specific heat values
    X = np.zeros(N_t)  # Susceptibility values
    Temperature = 0.4  # Fixed temperature to visualize the evolution of system

    # Boltzman constant is set to 1

    #   Algorithmic constants
    eq_step = 1024  # Number of Monte Carlo sweeps for equilibrium
    cal_step = 1024  # Number of Monte Carlo sweeps for calculation
    n1 = 1.0 / (cal_step * (N ** 2))
    n2 = 1.0 / ((cal_step ** 2) * (N ** 2))

    ###########################################################################################

    def __init__(self, lattice):
        '''
            Initialising class
        :param lattice: lattice dimension
        '''
        self.N = lattice

    def set_lattice(self, lattice):
        self.N = lattice

    def set_temperature(self, temperature):
        self.Temperature = temperature

    def set_temperature_point(self, n):
        self.N_t = n
        self.T = np.linspace(1.53, 3.28, self.N_t)
        self.E = np.zeros(self.N_t)  # Energy values
        self.M = np.zeros(self.N_t)  # Magnetization values
        self.C = np.zeros(self.N_t)  # Specific heat values
        self.X = np.zeros(self.N_t)  # Susceptibility values

    def set_temperature_range(self, temp1, temp2):
        if temp1 > temp2:
            m = temp1
            temp1 = temp2
            temp2 = m
        elif temp1 == temp2:
            temp2 += temp1

        self.T = np.linspace(temp1, temp2, self.N_t)

    def _initialize_state(self):
        '''
            Generates a random spin configuration for initial condition
        :param N: lattice dimension
        :return: 2-D array with spin configuration
        '''
        state = 2 * np.random.randint(2, size=(self.N, self.N)) - 1
        return state

    def monte_carlo(self, config, beta):
        '''
            Monte carlo algorithm

        :param config: current configuration of system
        :param beta: beta is defined by the inverse of the
                        product of Boltzman constant and temperature.
                        Boltzman constant is set to unity.
        :return: Accepted configuration after one monte carlo step
        '''
        N = self.N
        for i in range(N):
            for j in range(N):
                a = np.random.randint(0, N)
                b = np.random.randint(0, N)

                #   Randomly selecting a point on lattice
                s = config[a, b]

                #   Sum of neighboring points spin
                nb = config[(a + 1) % N, b] + config[a, (b + 1) % N] + config[(a - 1) % N, b] + config[a, (b - 1) % N]

                #   Delta Hamiltonian
                D_H = 2 * s * nb

                #   Configuration acceptance
                if D_H < 0:
                    s *= -1
                elif rand() < np.exp(-D_H * beta):  # rand() is between [0,1]
                    s *= -1  # rand() is w
                config[a, b] = s
        return config

    def calc_energy(self, config):
        '''
            Calculate energy of a given configuration
            Sum of {ij) S_i * S_j
        :param config: Configuration of system
        :return: Energy
        '''
        N = self.N
        energy = 0
        for i in range(len(config)):
            for j in range(len(config)):
                S = config[i, j]
                nb = config[(i + 1) % N, j] + config[i, (j + 1) % N] + config[(i - 1) % N, j] + config[i, (j - 1) % N]
                energy += -nb * S
        return energy / 4.

    def calc_magetization(self, config):
        '''
            Calculate normalized magnetization of a given configuration
            Sum of states of all points
        :param config: Configuration of system
        :return: Normalized magnetization
        '''
        mag = np.sum(config)
        return mag

    def calc_energy_T(self, temperature):
        '''
            Calculate expected value of energy in a specific temperature
        :param temperature: Given Temperature
        :return: <Energy> and <Energy^2>
        '''

        Energy = Energy2 = 0
        config = self._initialize_state()
        iT = 1.0 / temperature

        #   Equilibrium process
        for i in range(self.eq_step):
            self.monte_carlo(config, iT)  # Monte Carlo moves

        for i in range(self.cal_step):
            self.monte_carlo(config, iT)
            E = self.calc_energy(config)  # Calculate energy

            Energy += E
            Energy2 += E ** 2

        return Energy, Energy2

    def calc_magnet_T(self, temperature):
        '''
            Calculate expected value of magnetisation in a specific temperature
        :param temperature: Given Temperature
        :return: <M> and <M^2>
        '''

        Magnet = Magnet2 = 0
        config = self._initialize_state()
        iT = 1.0 / temperature
        iT2 = iT ** 2

        #   Equilibrium process
        for i in range(self.eq_step):
            self.monte_carlo(config, iT)  # Monte Carlo moves

        for i in range(self.cal_step):
            self.monte_carlo(config, iT)
            M = self.calc_magetization(config)  # Calculate the magnetisation

            Magnet += M
            Magnet2 += M ** 2

        return Magnet, Magnet2

    def calc_heat(self, temperature):
        '''
            Calculate specific heat in a given temperature
        :param temperature: Given temperature
        :return: <C>
        '''

        iT = 1.0 / temperature
        iT2 = iT ** 2

        Energy, Energy2 = self.calc_energy_T(temperature)

        return (self.n1 * Energy2 - self.n2 * Energy ** 2) * iT2

    def calc_suscep(self, temperature):
        '''
            Calculate susceptibility in a given temperature
        :param temperature: Given temperature
        :return: <X>
        '''

        iT = 1.0 / temperature

        Magnet, Magnet2 = self.calc_magnet_T(temperature)

        return (self.n1 * Magnet2 - self.n2 * Magnet ** 2) * iT

    def calc_EMCX(self):
        '''
            Calculate energy, magnetisation, specific heat and susceptibility
            with different temperatures
        :return:
        '''

        for index_t in range(self.N_t):
            iT = 1.0 / self.T[index_t]
            iT2 = iT ** 2
            Energy, Energy2 = self.calc_energy_T(self.T[index_t])
            Magnet, Magnet2 = self.calc_magnet_T(self.T[index_t])

            print("\nStep {} with Temperature {}\n".format(index_t, self.T[index_t]))

            self.E[index_t] = self.n1 * Energy
            print("Energy in temperature {} is {}".format(self.T[index_t], self.E[index_t]))

            self.M[index_t] = self.n1 * Magnet
            print("Magnetisation in temperature {} is {}".format(self.T[index_t], self.M[index_t]))

            self.C[index_t] = (self.n1 * Energy2 - self.n2 * Energy ** 2) * iT2
            print("Specific Heat in temperature {} is {}".format(self.T[index_t], self.C[index_t]))

            self.X[index_t] = (self.n1 * Magnet2 - self.n2 * Magnet ** 2) * iT
            print("Susceptibility in temperature {} is {}".format(self.T[index_t], self.X[index_t]))

    def plot_every(self):
        star = mpath.Path.unit_regular_star(6)
        circle = mpath.Path.unit_circle()
        # concatenate the circle with an internal cutout of the star
        verts = np.concatenate([circle.vertices, star.vertices[::-1, ...]])
        codes = np.concatenate([circle.codes, star.codes])
        cut_star = mpath.Path(verts, codes)

        ax = self.f.add_subplot(2, 2, 1)
        plt.scatter(self.T, self.E, s=10, marker='x', color='IndianRed')
        plt.xlabel("Temperature (T)", fontsize=20)
        plt.ylabel("Energy ", fontsize=20)
        plt.axis('tight')

        ax1 = self.f.add_subplot(2, 2, 2)
        plt.scatter(self.T, abs(self.M), s=10, marker='x', color='RoyalBlue')
        plt.xlabel("Temperature (T)", fontsize=20)
        plt.ylabel("Magnetisation ", fontsize=20)
        plt.axis('tight')

        ax2 = self.f.add_subplot(2, 2, 3)
        plt.scatter(self.T, self.C, s=10, marker='x', color='RoyalBlue')
        plt.xlabel("Temperature (T)", fontsize=20)
        plt.ylabel("Specific Heat ", fontsize=20)
        plt.axis('tight')

        ax3 = self.f.add_subplot(2, 2, 4)
        plt.scatter(self.T, self.X, s=10, marker='x', color='RoyalBlue')
        plt.xlabel("Temperature (T)", fontsize=20)
        plt.ylabel("Susceptibility ", fontsize=20)
        plt.axis('tight')

        plt.show()

    def plot_energy(self):
        '''
            Plot energy vs. temperature
        :return:
        '''
        # self.calc_EMCX()
        plt.scatter(self.T, self.E, s=20, marker='o', color='IndianRed')
        plt.xlabel("Temperature (T)", fontsize=20)
        plt.ylabel("Energy ", fontsize=20)
        plt.axis('tight')
        plt.show()

    def plot_magnetisation(self):
        '''
            Plot magnetisation vs. temperature
        :return:
        '''
        # self.calc_EMCX()
        plt.scatter(self.T, abs(self.M), s=20, marker='o', color='RoyalBlue')
        plt.xlabel("Temperature (T)", fontsize=20)
        plt.ylabel("Magnetisation ", fontsize=20)
        plt.axis('tight')
        plt.show()

    def plot_specific_heat(self):
        '''
            Plot specific heat vs. temperature
        :return:
        '''
        # self.calc_EMCX()
        plt.scatter(self.T, self.C, s=20, marker='o', color='RoyalBlue')
        plt.xlabel("Temperature (T)", fontsize=20)
        plt.ylabel("Specific Heat ", fontsize=20)
        plt.axis('tight')
        plt.show()

    def plot_susceptibility(self):
        '''
            Plot susceptibility vs. temperature
        :return:
        '''
        # self.calc_EMCX()
        plt.scatter(self.T, self.X, s=20, marker='o', color='RoyalBlue')
        plt.xlabel("Temperature (T)", fontsize=20)
        plt.ylabel("Susceptibility ", fontsize=20)
        plt.axis('tight')
        plt.show()

    def visualize(self, index, config, temp):
        N = self.N
        self.monte_carlo(config, 1.0 / temp)
        X, Y = np.meshgrid(range(N), range(N))
        plt.setp(self.ax.get_yticklabels(), visible=False)
        plt.setp(self.ax.get_xticklabels(), visible=False)
        self.ax.pcolormesh(X, Y, config, cmap=plt.get_cmap('RdBu'))
        plt.title('Time={}'.format(index))
        plt.axis('tight')

    def simulate(self, time, temperature=Temperature):
        '''
            Simulate the evolution of monte carlo algorithm
        :param temperature:
        :return:
        '''

        config = self._initialize_state()

        ani = anim.FuncAnimation(self.f, self.visualize, time,
                                 fargs=(config, temperature),
                                 interval=1, blit=False)
        plt.show()


if __name__ == '__main__':

    '''
        Ising model with 16x16 lattice and 50 temperature points to  plot
        Simulation time is et to 1000
    '''
    ising = Ising(16)

    eq = int(input("Simulation or plot \n "
                   "\t[1] Simulation\n"
                   "\t[2] Plot\n"
                   "\tSelect (1/2): "))
    if eq == 2:
        temp_poinr = int(input("\nTemperature Point: "))
        ising.set_temperature_point(temp_poinr)
        ising.calc_EMCX()
        ising.plot_every()
    elif eq == 1:
        temp = float(input("\nTemperature : "))
        if temp < 0: temp = -temp
        ising.simulate(1000, temp)
