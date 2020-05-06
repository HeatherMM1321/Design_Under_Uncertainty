# Heather Miller
# ME615 HW3
# started 5/6/20

import numpy as np
from scipy.stats import norm

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

def calculate_weight(diameter, length, rho_cu, rho_fe):
    # calculate the weight of the motor
    weight = (rho_cu * (Awa*Lwa + Awf*Lwf) + rho_fe * length * np.pi * pow(diameter + ds, 2)) - 22
    return weight

# variable info
variables = {"diameter": (0.075, 0.005),
             "length": (0.095, 0.005),
             "rho_cu": (8.94e3, 100),
             "rho_fe": (7.98e3, 100)}


def finite_diff(variables):
    #compute derivatives of funtion from finite difference approximation

    # These are the mu and sigmas of each variable
    D = (variables["diameter"][0], variables["diameter"][1])
    L = (variables["length"][0], variables["length"][1])
    rho_cu = (variables["rho_cu"][0], variables["rho_cu"][1])
    rho_fe = (variables["rho_fe"][0], variables["rho_fe"][1])

    d_dD = (calculate_weight(D[0] + 0.1 * D[1], L[0], rho_cu[0], rho_fe[0]) -
            calculate_weight(D[0] - 0.1 * D[1], L[0], rho_cu[0], rho_fe[0])) / (0.2 * D[1])

    d_dL = (calculate_weight(D[0], L[0] + 0.1 * L[1], rho_cu[0], rho_fe[0]) -
            calculate_weight(D[0], L[0] - 0.1 * L[1], rho_cu[0], rho_fe[0])) / (0.2 * L[1])

    d_dcu = (calculate_weight(D[0], L[0], rho_cu[0] + 0.1 * rho_cu[1], rho_fe[0]) -
             calculate_weight(D[0], L[0], rho_cu[0] - 0.1 * rho_cu[1], rho_fe[0])) / (0.2 * rho_cu[1])

    d_dfe = (calculate_weight(D[0], L[0], rho_cu[0], rho_fe[0] + 0.1 * rho_fe[1]) -
             calculate_weight(D[0], L[0], rho_cu[0], rho_fe[0] - 0.1 * rho_fe[1])) / (0.2 * rho_fe[1])

    return [d_dD, d_dL, d_dcu, d_dfe]

def update_eq(u, variables):

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


if __name__ == "__main__":
    u = [0, 0, 0, 0]
    diff = 1
    k = 1
    while diff > 0.001 and k < 20:
        uk = update_eq(u, variables)
        diff = abs(np.linalg.norm(uk) - np.linalg.norm(u))
        u = uk
        prob = norm.cdf(np.linalg.norm(u))

        print("Iteration #", k, ":\n", u, "\nIntermediate Probability", prob, "\nIntermediate Beta",
              np.linalg.norm(u), "\n")
        k += 1

