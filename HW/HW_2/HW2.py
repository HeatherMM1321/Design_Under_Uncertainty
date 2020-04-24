# Heather Miller
# ME615 HW2
# started 4/22/20

import numpy as np
from scipy.stats import norm


# Motor Design problem : weight of motor must be <22kgs

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
ap = p*1  # parallel circuit paths (equals poles)
eff = 0.85  # efficiency
bfc = 40e-3  # pole depth (m)
fcf = 0.55  # field coil spacing factor
Awa = 2.0074e-006  # cross sectional area of armature winding (m^2)
Awf = 0.2749e-006  # cross sectional area of field coil winding (m^2)

#########################################################################################

samples = 100000  # number of samples

def calculate_weight(diameter, length, rho_cu, rho_fe):
    # calculate the weight of the engine
    weight = rho_cu * (Awa*Lwa + Awf*Lwf) + rho_fe * length * np.pi * pow(diameter + ds, 2)
    return weight

def prob_success_mc(weights, limit):
    # determine the probability of success that the weight of engine will be less than a limit
    # add 1 to limit_sum every time weight is under limit weight
    x = len(weights)
    limit_sum = 0
    for i in weights:
        if i < limit:
            limit_sum += 1
    return limit_sum/x

########################################################################################
# Part 1 design parameters: Using Monte Carlo and Normal Distribution
# D~N(7.5, 0.5) %m rotor diameter
# L~N(9.5, 0.5) %m rotor axial length
# dcu~N(8.94e3, 100) %copper density density at 25C (kg/m^3)
# dfe~N(7.98e3, 100) %iron density density at 25C (kg/m^3)

d_n = np.random.normal(0.075, 0.05, samples)  # rotor diameter (cm)
l_n = np.random.normal(0.095, 0.05, samples)  # rotor axial length (cm)
rho_cu_n = np.random.normal(8.94e3, 100, samples)  # copper density density at 25C (kg/m^3)
rho_fe_n = np.random.normal(7.98e3, 100, samples)  # iron density density at 25C (kg/m^3)
weight_p1 = [calculate_weight(d_n[i], l_n[i], rho_cu_n[i], rho_fe_n[i]) for i in range(samples)]
prob_success = prob_success_mc(weight_p1, 22)
print("The probability of a motor being less than 22kg with the parameters in #1 is", prob_success, "%")


#######################################################################################################
# Part 2 design parameters: Using Monte Carlo and Uniform Distribution
# D ~ Uniform(6.5, 8.5)
# L ~ Uniform (8.5, 10.5)
# dcu ~ Uniform (8840, 9040)
# dfe ~ Uniform (7880, 8080)
d_u = np.random.uniform(0.065, 0.085, samples)  # rotor diameter (cm)
l_u = np.random.uniform(0.085, 0.105, samples)  # rotor axial length (cm)
rho_cu_u = np.random.uniform(8840, 9040, samples)  # copper density density at 25C (kg/m^3)
rho_fe_u = np.random.uniform(7880, 8080, samples)  # iron density density at 25C (kg/m^3)
weight_p2 = [calculate_weight(d_u[i], l_u[i], rho_cu_u[i], rho_fe_u[i]) for i in range(samples)]
prob_success = prob_success_mc(weight_p2, 22)
print("The probability of a motor being less than 22kg with the parameters in #2 is", prob_success, "%")


###########################################################################################################
# Part 3 design parameters: Using MVFOSM method##
#(mu, sigma, first derivative)
diameter = (.075, 0.5) # m rotor diameter
length = (0.095, 0.05) # m rotor axial length
rho_cu = (8.94e3, 100) # copper density density at 25C (kg/m^3)
rho_fe = (7.98e3, 100) # iron density density at 25C (kg/m^3)

function_mu = calculate_weight(diameter[0], length[0], rho_cu[0], rho_fe[0])

diameter_output = (calculate_weight(diameter[0]+(0.1 * diameter[1]), length[0], rho_cu[0], rho_fe[0]) -
         calculate_weight(diameter[0]-(0.1 * diameter[1]), length[0], rho_cu[0], rho_fe[0]))/ 2*(0.1 * diameter[1])

length_output = (calculate_weight(diameter[0], length[0]+(0.1 * length[1]), rho_cu[0], rho_fe[0]) - \
         calculate_weight(diameter[0], length[0]-(0.1 * length[1]), rho_cu[0], rho_fe[0]))/ 2*(0.1 * length[1])

cu_output = (calculate_weight(diameter[0], length[0], rho_cu[0]+(0.1 * rho_cu[1]), rho_fe[0]) - \
         calculate_weight(diameter[0], length[0], rho_cu[0]-(0.1 * rho_cu[1]), rho_fe[0]))/ 2*(0.1 * rho_cu[1])

fe_output = (calculate_weight(diameter[0], length[0], rho_cu[0], rho_fe[0]+(0.1 * rho_fe[1])) - \
         calculate_weight(diameter[0], length[0], rho_cu[0], rho_fe[0]-(0.1 * rho_fe[1])))/ 2*(0.1 * rho_fe[1])


variable_derivatives = [diameter_output, length_output, cu_output, fe_output]
variable_sigmas = [diameter[1], length[1], rho_cu[1], rho_fe[1]]


#no correlation
correlation_matrix = [[1, 0, 0 ,0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]

def calculate_sigma(variable_sigmas, variable_derivatives, correlation_matrix):
    sigma_squared = 0
    for i in  range(len(variable_derivatives)):
        for j in range(len(variable_derivatives)):
            sigma_squared += variable_derivatives[i] * variable_derivatives[j] * correlation_matrix[i][j] \
                             * variable_sigmas[i] * variable_sigmas[j]
    return np.sqrt(sigma_squared)

sigma = calculate_sigma(variable_sigmas, variable_derivatives, correlation_matrix)

success_probability = norm.cdf(22, function_mu, sigma)
print('The probability of a motor being less than 22kg with the parameters in #3 is', success_probability, '%')


########################################################################################################
# Part 4: Repeat part 3 but with correlation matirx##
# how does the correlation change the solution??

correlation_matrix = [[1, .2, .3, .7], [.2, 1, .5, .6], [.3, .5, 1, .2], [.7, .6, .2, 1]]
sigma = calculate_sigma(variable_sigmas, variable_derivatives, correlation_matrix)

success_probability = norm.cdf(22, function_mu, sigma)
print('The probability of a motor being less than 22kg with the parameters in #4 is', success_probability, '%')
