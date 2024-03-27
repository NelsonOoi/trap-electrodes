import numpy as np
import scipy.optimize as opt
import scipy.constants as ct
import matplotlib.pyplot as plt

'''
#############################################

        Inspired by Luke Qi's thesis.

#############################################
'''

# ion_charges = np.array([1., 1., 1., 1., 1.]) * ct.elementary_charge
ions_charges = np.ones(9) * ct.elementary_charge
# ion_charges[1] = 1 * ct.elementary_charge
# ion_charges = np.array([1.])
# ion_charges = np.array([1., 1.]) * ct.elementary_charge
n_ions = len(ions_charges)
print(n_ions)
# target_interion_spacing = 5.
target_interion_spacing = 5e-6

target_interion_spacings_arr = np.ones(n_ions - 1) * target_interion_spacing
position_ansatz = np.arange(n_ions) * target_interion_spacing
recentering = (n_ions - 1) * target_interion_spacing
position_ansatz = position_ansatz - recentering / 2

global cascaded_position_ansatz
cascaded_position_ansatz = position_ansatz
ab_ansatz = [1e6, 1e16]
ab_bounds = ((1e4, 1e20), (1e2, 1e24))
# ab_ansatz = [1e-12, 1e-32]
# ab_bounds = ((0, 1e-8), (0, 1e-8))
# ab_ansatz = [1e-6, 1e-8]
# ab_bounds = ((0, 1e-4), (0, 1e-4))

# pos_bounds = ((None, 0), (None, 0), (0, None), (0, None))
# pos_bounds = ((None, 0), (0, None))

print(position_ansatz)

# Cost function for outer optimization loop

def eq_positions(ab, cascaded_position_ansatz=cascaded_position_ansatz):
    # ab = [2.588e-6, 2.346e-8]
    # ct.epsilon_0 has units A2·s4·kg−1·m−3 -> 1e-18 • A2·s4·kg−1·µm−3
    # Cost function for inner optimization loop 
    # e0_scaling = 1e-18
    # print(k)
    # scaling = 1e-6
    print(ab)
    def scaled_potential(ions_pos):
        energy_scaling = 1e19 #prevents the x^4 term from vanishing
        e0_scaling = 1.
        n_ions = len(ions_pos)
        k = 1 / (4 * np.pi * ct.epsilon_0 * e0_scaling) * 1/2
        u = 0
        for i in range(n_ions):
            u += ions_charges[i] * (ab[0] * (ions_pos[i]) **2 + ab[1] * (ions_pos[i]) **4) * energy_scaling
            interaction_potential = 0
            other_ions = list(range(n_ions))
            other_ions.remove(i)
            for j in other_ions:
                interaction_potential += ions_charges[i] * k * ions_charges[j] / (np.abs(ions_pos[j] - ions_pos[i])) * energy_scaling
            u += interaction_potential
        return u

    # inner_opt_options = {
    #                     #  'disp':True,
    #                      'maxiter': int(1e3),
    #                      'maxfev': int(1e3),
    #                     #  'return_all': True,
    #                     #  'tol': 1e-17, #cobyla
    #                     #  'rhobeg': 0.1, #cobyla
    #                     #  'ftol': 0.00001,
    #                     #  'fatol': 1e-10,
    #                     #  'xatol': 1e-10,
    #                     #  'xrtol': 0.001,
    #                     #  'norm': 2,
    #                     'eps': 0.0001, #bfgs
    #                     # 'gtol': 1e-7, #l-bfgs-b
    #                     'adaptive': True,
    #                      }
    # # inner_opt_method = 'Nelder-Mead'
    # inner_opt_method = 'BFGS'
    # # inner_opt_method = 'L-BFGS-B'
    # # inner_opt_method = 'CG'
    # # inner_opt_method = 'COBYLA'
    # # inner_opt_method = 'dogleg'
    # res_potential = opt.minimize(fun=scaled_potential,
    #                              x0=cascaded_position_ansatz,
    #                             method=inner_opt_method,
    #                             options=inner_opt_options,
    #                         #  bounds=pos_bounds
    #                          )
    # minimizer_kwargs = {'method': 'Nelder-Mead'}
    minimizer_kwargs = {'method': 'BFGS'}
    res_potential = opt.basinhopping(func=scaled_potential,
                                 x0=cascaded_position_ansatz,
                                 niter=int(10),
                                 minimizer_kwargs=minimizer_kwargs,
                                #  disp=True,
                                #  method=inner_opt_method,
                                #  options=inner_opt_options,
                            #  bounds=pos_bounds
                                )
    # res_potential = opt.basinhopping(func=potential, x0=position_ansatz)
    print(res_potential)
    positions = res_potential.x
    positions = np.sort(positions)
    interion_spacings = np.diff(positions)
    cascaded_position_ansatz = positions
    # takes difference between each successive position
    # print(np.abs(interion_spacings) - target_interion_spacings_arr)
    # print((np.abs(interion_spacings) - target_interion_spacings_arr)**2)

    # cost = np.sum((np.abs(interion_spacings) - target_interion_spacings_arr)**2)
    # cost = cost * 1e12
    # print(positions, interion_spacings, cost)
    # return cost
    return positions, interion_spacings

def squared_dist_deviation(ab):
    positions, interion_spacings = eq_positions(ab)
    cost = np.sum((np.abs(interion_spacings) - target_interion_spacings_arr)**2)
    cost = cost * 1e6
    # print(positions, interion_spacings, cost)
    return cost


outer_opt_options = {
                     'disp':True,
                     'maxiter': int(1e2),
                    #  'maxfev': int(1e8),
                    #  'tol': 0.00001,
                    #  'gtol': 1e-7,
                     'fatol': 1e-10,
                     'xatol': 1e5,
                     'adaptive': True,
                    #  'xrtol': 0.001
                    'eps': 1e4,
                     }
outer_opt_method = 'Nelder-Mead'
# outer_opt_method = 'L-BFGS-B'
# outer_opt_method = 'BFGS'
res = opt.minimize(fun=squared_dist_deviation,
                    x0=ab_ansatz,
                    method=outer_opt_method,
                    options=outer_opt_options,
                    bounds=ab_bounds,
               )
# minimizer_kwargs = {'method': 'BFGS'}
# res = opt.basinhopping(func=squared_dist_deviation,
#                     x0=ab_ansatz,
#                     niter=10,
#                     disp=True,
#                     # method=outer_opt_method,
#                     # options=outer_opt_options,
#                     # bounds=ab_bounds,
#                     minimizer_kwargs=minimizer_kwargs,
#                )
# res = opt.shgo(squared_dist_deviation, ab_bounds)
# res = opt.basinhopping(squared_dist_deviation, x0=ab_ansatz)

# Outer optimization loop
    
# res = minimize(squared_dist_deviation, x0, method=method, options=options, bounds=bounds)

# print(squared_dist_deviation([1e4, 1]))
# print(squared_dist_deviation([1e-5, 2e-10]))
# print(squared_dist_deviation([2.588e-6, 2.346e-8]))
# print(squared_dist_deviation([2e6, 0]))
# print(squared_dist_deviation([2.588e12, 0]))
# print(squared_dist_deviation([2.588e-18, 2.346e-32]))
# print(squared_dist_deviation([2e-17, 2.346e-20]))

def potential(ab, x, ions_charges, ions_pos, current_ion = 0):
    # energy_scaling = 1e12 #prevents the x^4 term from vanishing
    # e0_scaling = 1.
    k = 1 / (4 * np.pi * ct.epsilon_0) * 1/2
    u = 0
    n_ions = len(ions_charges)
    u += (ab[0] * (x) **2 + ab[1] * (x) **4)
    interaction_potential = 0
    other_ions = list(range(n_ions))
    other_ions.remove(current_ion)
    for j in other_ions:
        interaction_potential += k * ions_charges[j] / (np.abs(ions_pos[j] - x))
    u += interaction_potential
    return u

mu = 1e-6

res_ab = res.x
plt.figure(figsize=(10,6))
res_pos, res_spacings = eq_positions(res_ab)
ions_potentials = []
print(res_pos, res_spacings, res_ab)
for ab in res_ab:
    print(ab)

for i in range(n_ions):
    delta = target_interion_spacing * 0.75
    x_around_ion = np.linspace(res_pos[i] - delta, res_pos[i] + delta)
    potential_around_ion = potential(res_ab, x_around_ion, ions_charges, res_pos, current_ion=i)
    plt.plot(x_around_ion / mu, potential_around_ion)
    print(np.min(potential_around_ion))
    ions_potentials.append(potential(res_ab, res_pos[i], ions_charges, res_pos, current_ion=i))
print(ions_potentials)

plt.scatter(x=res_pos/mu, y=ions_potentials, color='r', label=f'$40Ca^+$ ions')

x = np.linspace(res_pos[0]-5e-7, res_pos[-1]+5e-7)
curve = res_ab[0] * x**2 + res_ab[1] * x**4
curve_mean = np.mean(curve)
plt.plot(x / mu, curve, label='Trap potential')
plt.ylabel('V (V)')
plt.xlabel('x (µm)')
y = np.ones(len(res_pos)) * curve_mean
# plt.scatter(x=res_pos/mu, y=y, color='r', label=f'$40Ca^+$ ions')
plt.legend()
plt.show()
