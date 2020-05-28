# Heather Miller
# ME615 HW3
# started 5/6/20

import numpy as np
from scipy.stats import norm
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss
from pyDOE2 import fullfact
from pears_cdf import *

####################################################
# Given Constants
# Motor Design Variables
Lwa = 14.12  # armature wire length (m)
Lwf = 309.45  # field wire length (m)
ds = 0.00612  # slot depth (m)

# Coupling Variables b, shared with control problem
n = 122  # rotational speed (rev/s)
v = 40  # design voltage (V)
pmin = 3.94  # minimum required power (kW)
ymin = 5.12e-3  # minimum required torque (kNm)

# Parameter Vector a (constants)
fi = 0.7  # pole arc to pole pitch ratio
p = 2  # number of poles
s = 27  # number of slots (teeth on rotor)
rho = 1.8e-8  # resistivity (ohm-m) of copper at 25C

# Derived parameters and constants
mu0 = 4 * np.pi * 1e-7  # magnetic constant
ap = p  # parallel circuit paths (equals poles)
eff = 0.85  # efficiency
bfc = 40e-3  # pole depth (m)
fcf = 0.55  # field coil spacing factor
Awa = 2.0074e-006  # cross sectional area of armature winding (m^2)
Awf = 0.2749e-6  # cross sectional area of field coil winding (m^2)

#########################################################################################
# Part 1  : Using FORM


def calculate_weight(diameter, length, rho_cu, rho_fe, limit=22):
    # calculate the weight of the motor
    # set limit=0 you are not calculating in u space
    weight = (rho_cu * (Awa*Lwa + Awf*Lwf) + rho_fe * length * np.pi * pow(diameter + ds, 2)) - limit
    return weight


def finite_diff(variables):
    # returns the first derivative of diameter, length, rho_cu, and rho_fe

    # These are the mu and sigmas of each variable
    D = (variables["diameter"][0], variables["diameter"][1])
    L = (variables["length"][0], variables["length"][1])
    rho_cu = (variables["rho_cu"][0], variables["rho_cu"][1])
    rho_fe = (variables["rho_fe"][0], variables["rho_fe"][1])

    # first derivative of each variable
    d_dD = (calculate_weight(D[0] + 0.1 * D[1], L[0], rho_cu[0], rho_fe[0]) -
            calculate_weight(D[0] - 0.1 * D[1], L[0], rho_cu[0], rho_fe[0])) / (0.2 * D[1])

    d_dL = (calculate_weight(D[0], L[0] + 0.1 * L[1], rho_cu[0], rho_fe[0]) -
            calculate_weight(D[0], L[0] - 0.1 * L[1], rho_cu[0], rho_fe[0])) / (0.2 * L[1])

    d_dcu = (calculate_weight(D[0], L[0], rho_cu[0] + 0.1 * rho_cu[1], rho_fe[0]) -
             calculate_weight(D[0], L[0], rho_cu[0] - 0.1 * rho_cu[1], rho_fe[0])) / (0.2 * rho_cu[1])

    d_dfe = (calculate_weight(D[0], L[0], rho_cu[0], rho_fe[0] + 0.1 * rho_fe[1]) -
             calculate_weight(D[0], L[0], rho_cu[0], rho_fe[0] - 0.1 * rho_fe[1])) / (0.2 * rho_fe[1])

    return [d_dD, d_dL, d_dcu, d_dfe]


def u_update(u, variables):
    # returns the updated u

    variable_mus = [variables["diameter"][0], variables["length"][0], variables["rho_cu"][0],
                               variables["rho_fe"][0]]
    variable_sigmas = [variables["diameter"][1], variables["length"][1], variables["rho_cu"][1],
                               variables["rho_fe"][1]]

    # this is iterating through each of the variables mu and sigma and converting it into u space
    x = [(variable_mus[i] + u[i]*variable_sigmas[i]) for i in range(len(u))]

    # this is iterating through the first derivatives of each of the variables and multiplying them by their respective
    # sigmas
    grad_g = [(finite_diff(variables)[i]* variable_sigmas[i]) for i in range(len(u))]

    # this is plugging in the u space version of each variable into the equation
    g_u = calculate_weight(x[0], x[1], x[2], x[3])

    # updating u
    u_k = [((np.dot(grad_g, u) - g_u)/(np.dot(grad_g, grad_g))) * i for i in grad_g]

    return u_k


def find_reliability(u, diff_threshold, iterations):
    # search for the reliability of function
    # returns the number of iterations until convergence, the final beta, and the probability of success

    diff = 1
    k = 1
    while diff > diff_threshold and k < iterations:
        uk = u_update(u, variables)
        diff = abs(np.linalg.norm(uk) - np.linalg.norm(u))
        u = uk
        beta = np.linalg.norm(u)
        prob = norm.cdf(np.linalg.norm(beta))

        k += 1
    return k-1, beta, prob


#####################################################################################################
# Part 2 : Using Full Tensor Numerical Integration Method
# 3 nodes per uncertain input
# same uncertainties as problem 1

def convert_norm(number_of_nodes, variables):
    # convert Gauss–Hermite nodes and weights to integrate specific distributions

    variable_mus = [variables["diameter"][0], variables["length"][0], variables["rho_cu"][0],
                               variables["rho_fe"][0]]
    variable_sigmas = [variables["diameter"][1], variables["length"][1], variables["rho_cu"][1],
                               variables["rho_fe"][1]]

    points, weights = hermgauss(number_of_nodes)
    variable_points =[]
    variable_weights = [i / np.sqrt(np.pi) for i in weights]

    for j in range(len(variables)):
        variable_points.append([(np.sqrt(2) * i * variable_sigmas[j]) + variable_mus[j] for i in points])

    return variable_points, variable_weights


def convert_uniform(number_of_nodes, variables):
    # convert Gauss-Laguerre nodes and weights to integrate specific distributions

    variable_lower = [variables["diameter"][0], variables["length"][0], variables["rho_cu"][0],
                               variables["rho_fe"][0]]
    variable_upper = [variables["diameter"][1], variables["length"][1], variables["rho_cu"][1],
                               variables["rho_fe"][1]]

    points, weights = leggauss(number_of_nodes)
    variable_points =[]
    variable_weights = [i / 2 for i in weights]

    for j in range(len(variables)):
        variable_points.append([(variable_upper[j] - variable_lower[j]) * i / 2 +
                                (variable_upper[j] + variable_lower[j]) / 2 for i in points])

    return variable_points, variable_weights


def full_factorial(points, weights):
    # perform a full factorial on the converted points and weights

    D = points[0]
    L = points[1]
    cu = points[2]
    fe = points[3]

    g_points = []
    big_W = []

    indx = fullfact([3, 3, 3, 3])

    for i in indx:
        idx_D = int(i[0])
        idx_L = int(i[1])
        idx_cu = int(i[2])
        idx_fe = int(i[3])
        g_points.append(calculate_weight(D[idx_D], L[idx_L], cu[idx_cu], fe[idx_fe], limit=0))
        big_W.append(weights[idx_D]*weights[idx_L]*weights[idx_cu]*weights[idx_fe])
        # table_info.append([D[idx_D], L[idx_L], cu[idx_cu], fe[idx_fe], g_point, weight])

    # df_outputs = pd.DataFrame(table_info, columns=['D', 'L', 'cu', 'fe', 'G(pts)', 'W'])
    return g_points, big_W


def cal_sys_moments(g_points, big_W):
    # calculate the moments of the function

    g_mean = np.dot(g_points, big_W)
    g_sigma = np.sqrt(np.dot((g_points-g_mean)**2, big_W))
    g_skew = (np.dot((((g_points-g_mean)**3)/g_sigma**3), big_W))**2
    g_kurtosis = np.dot((((g_points-g_mean)**4)/g_sigma**4), big_W)

    return g_mean, g_sigma, g_skew, g_kurtosis


if __name__ == "__main__":
    # Compute all HW answers

    # Part 1: Normal Distribution (FORM)
    variables = {"diameter": (0.075, 0.005),
                 "length": (0.095, 0.005),
                 "rho_cu": (8.94e3, 100),
                 "rho_fe": (7.98e3, 100)}
    u = [0, 0, 0, 0]
    k, beta, prob = find_reliability(u, .001, 20)

    print("Problem 1:\n# of Iterations:", k, "\nBeta:", beta, "\nProbability of Success:", prob)

    # Part 2: Normal Distribution (Full Tensor Numerical Integration Method)
    points, weights = convert_norm(3, variables)
    g_points, big_W = full_factorial(points, weights)
    g_mean, g_sigma, g_skew, g_kurtosis = cal_sys_moments(g_points, big_W)
    p, type = pearson_fit([-np.inf, 22], g_mean, g_sigma, g_skew, g_kurtosis)

    print("\nProblem 2:\nMean:",g_mean,"\nStandard Deviation:", g_sigma, "\nSkew:", g_skew, "\nKurtosis:", g_kurtosis,
          "\nProbaility of Success:", p)

    # Part 3: Uniform Distribution (Full Tensor Numerical Integration Method)
    variables_3 = {"diameter": (0.065, 0.085),
                   "length": (0.085, 0.105),
                    "rho_cu": (8840, 9040),
                    "rho_fe": (7880, 8080)}

    points, weights = convert_uniform(3, variables_3)
    g_points, big_W = full_factorial(points, weights)
    g_mean, g_sigma, g_skew, g_kurtosis = cal_sys_moments(g_points, big_W)
    p, type = pearson_fit([-np.inf, 22], g_mean, g_sigma, g_skew, g_kurtosis)

    print("\nProblem 3:\nMean:",g_mean,"\nStandard Deviation:", g_sigma, "\nSkew:", g_skew, "\nKurtosis:", g_kurtosis,

          "\nProbaility of Success:", p)

