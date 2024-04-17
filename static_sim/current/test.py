from trapsim import *
from electrode import (System, PolygonPixelElectrode, euler_matrix,
    PointPixelElectrode, PotentialObjective,
    PatternRangeConstraint, shaped, utils)
import shapely.geometry as sg
import shapely.ops as so
import shapely as sh

s, electrodes, electrodes_dict = load_trap(filename='single_chip.gds',
            electrode_layer=37,
            ito_layer=12, electrode_mapping=default_electrode_mapping,
            electrode_ordering=electrode_ordering, plot=False,
            xlim=(-5000,5000), ylim=(-3000,3000), trap_center=default_trap_center, buildup=False)

dc_axial_set_file = 'Apr15_gds_Vs_axial.csv'
dc_tilt_set_file = 'Apr16_approx_Vs_tilt.csv'
dc_gds_tilt_set_file = 'Apr16_gds_Vs_tilt.csv'
dc_tilt_karan_file = 'Vs_tilt_karan.csv'

dc_axial_set, dc_tilt_set, dc_tilt_karan, dc_gds_tilt_set = read_electrode_voltages([dc_axial_set_file, dc_tilt_set_file, dc_tilt_karan_file, dc_gds_tilt_set_file])

o_traprf = 2*np.pi * 35e6
q = 1 * ct.elementary_charge
m = 40 * ct.atomic_mass
l = 1e-6
vp = 20.

with s.with_voltages(dc_tilt_karan):
    s["r"].rf = vp*np.sqrt(q/m)/(2*l*o_traprf)
    # plot_potential(s=s)
    plot_field(s=s, is_dc_potential=True, x_grid_bounds=(-10,10), y_grid_bounds=(-10,10), z_grid_bounds=(40,60))

with s.with_voltages(dc_tilt_set):
    s["r"].rf = vp*np.sqrt(q/m)/(2*l*o_traprf)
    # plot_potential(s=s)
    plot_field(s=s, is_dc_potential=True, x_grid_bounds=(-10,10), y_grid_bounds=(-10,10), z_grid_bounds=(40,60))

with s.with_voltages(dc_gds_tilt_set):
    s["r"].rf = vp*np.sqrt(q/m)/(2*l*o_traprf)
    # plot_potential(s=s)
    plot_field(s=s, is_dc_potential=True, x_grid_bounds=(-10,10), y_grid_bounds=(-10,10), z_grid_bounds=(40,60))

plt.show()
# shape = [[ -560. ,   102.5],
#        [ -560. ,  1112.5],
#        [-1135. ,  1595. ],
#        [-1135. ,  1745. ],
#        [ -985. ,  1745. ],
#        [ -985. ,  1585. ],
#        [ -440. ,  1102.5],
#        [ -440. ,   102.5],
#        [ -560. ,   102.5]]

# electrodes = [( '2',
#     [shape]
# )]


# s = System([PolygonPixelElectrode(name=n, paths=map(np.array, p))
#                 for n, p in electrodes])

# with s.with_voltages(dcs=[1.]):
#     print(s.electrical_potential([0., 0., 51.7], 'dc'))

# shape1 = sg.Polygon(shape)
# shape1 = sg.polygon.orient(shape1, sign=1.0)
# el_coordinates = np.array(sg.mapping(shape1).get('coordinates')[0])

# electrodes1 = [( '2',
#     [el_coordinates]
# )]

# print(shape1)

# s1 = System([PolygonPixelElectrode(name=n, paths=map(np.array, p))
#                 for n, p in electrodes1])

# with s1.with_voltages(dcs=[1.]):
#     print(s1.electrical_potential([0., 0., 51.7], 'dc'))
# files = ['Vs_axial_karan.csv', 'Vs_axial_karan2.csv', 'Vs_axial125_karan_recentered.csv', 'Vs_tilt_karan.csv']
# overall = read_electrode_voltages(files=files)

# for voltage_set in overall:
#     print(voltage_set)

# import numpy as np
# import scipy.optimize as opt
# import scipy.constants as ct
# import matplotlib.pyplot as plt

# ions_charges = np.ones(9) * ct.elementary_charge
# n_ions = len(ions_charges)
# # target_interion_spacing = 5.
# target_interion_spacing = 5e-6

# target_interion_spacings_arr = np.ones(n_ions - 1) * target_interion_spacing
# position_ansatz = np.arange(n_ions) * target_interion_spacing
# recentering = (n_ions - 1) * target_interion_spacing
# position_ansatz = position_ansatz - recentering / 2
# # ab_ansatz = [1e6, 1e16]
# # ab_bounds = ((0, 1e8), (0, 1e24))
# # ab_ansatz = [1e-12, 1e-32]
# # ab_bounds = ((0, 1e-8), (0, 1e-8))
# # ab_ansatz = [1e-6, 1e-8]
# ab_bounds = ((0, 1e10), (0, 1e20))
# # ab = [5.39518890e+05, 6.93869436e+15]
# # ab = [4.49953968e+05, 4.27092958e+15]
# ab = [4.51360317e+05, 4.16195171e+15]
# # ab = [449953.968, 4270929580000000.0]

# # ab = [2.588e-6, 2.346e-8]
# # print(ab)
# # Cost function for inner optimization loop 
# def scaled_potential(ions_pos):
#     energy_scaling = 1e19 #prevents the x^4 term from vanishing
#     e0_scaling = 1.
#     n_ions = len(ions_pos)
#     k = 1 / (4 * np.pi * ct.epsilon_0 * e0_scaling) * 1/2
#     u = 0
#     for i in range(n_ions):
#         u += ions_charges[i] * (ab[0] * (ions_pos[i]) **2 + ab[1] * (ions_pos[i]) **4) * energy_scaling
#         interaction_potential = 0
#         other_ions = list(range(n_ions))
#         other_ions.remove(i)
#         for j in other_ions:
#             interaction_potential += ions_charges[i] * k * ions_charges[j] / (np.abs(ions_pos[j] - ions_pos[i])) * energy_scaling
#         u += interaction_potential
#     return u

# inner_opt_options = {
#                     #  'disp':True,
#                         'maxiter': int(1e6),
#                         'maxfev': int(1e3),
#                     #  'return_all': True,
#                     #  'tol': 1e-17, #cobyla
#                     #  'rhobeg': 0.1, #cobyla
#                     #  'ftol': 0.00001,
#                         'fatol': 1e-10,
#                         'xatol': 1e-10,
#                     #  'xrtol': 0.001,
#                     #  'norm': 2,
#                     # 'eps': 0.0001, #bfgs
#                     # 'gtol': 1e-7, #l-bfgs-b
#                     'adaptive': True,
#                         }
# inner_opt_method = 'Nelder-Mead'
# # inner_opt_method = 'BFGS'
# # inner_opt_method = 'L-BFGS-B'
# # inner_opt_method = 'CG'
# # inner_opt_method = 'COBYLA'
# # inner_opt_method = 'dogleg'
# # res_potential = opt.minimize(fun=scaled_potential,
# #                              x0=position_ansatz,
# #                             method=inner_opt_method,
# #                             options=inner_opt_options,
# #                         #  bounds=pos_bounds
# #                          )
# # minimizer_kwargs = {'method': 'Nelder-Mead'}
# minimizer_kwargs = {'method': 'BFGS'}
# res_potential = opt.basinhopping(func=scaled_potential,
#                                 x0=position_ansatz,
#                                 niter=int(10),
#                                 disp=True,
#                                 minimizer_kwargs=minimizer_kwargs
# )

# def potential(ab, x, ions_charges, ions_pos, current_ion = 0):
#     # energy_scaling = 1e12 #prevents the x^4 term from vanishing
#     # e0_scaling = 1.
#     k = 1 / (4 * np.pi * ct.epsilon_0) * 1/2
#     u = 0
#     n_ions = len(ions_charges)
#     u += ions_charges[current_ion] * (ab[0] * (x) **2 + ab[1] * (x) **4)
#     interaction_potential = 0
#     other_ions = list(range(n_ions))
#     other_ions.remove(current_ion)
#     for j in other_ions:
#         interaction_potential += ions_charges[current_ion] * k * ions_charges[j] / (np.abs(ions_pos[j] - x))
#     u += interaction_potential
#     return u

# mu = 1e-6

# plt.figure(figsize=(10,6))
# res_ab = ab
# res_pos = np.sort(res_potential.x)
# res_spacings = np.diff(res_pos)
# ions_potentials = []
# print(res_pos, res_spacings, res_ab)

# for i in range(n_ions):
#     delta = np.min(res_spacings) * 0.5
#     x_around_ion = np.linspace(res_pos[i] - delta, res_pos[i] + delta)
#     potential_around_ion = potential(res_ab, x_around_ion, ions_charges, res_pos, current_ion=i)
#     plt.plot(x_around_ion / mu, potential_around_ion)
#     print(np.min(potential_around_ion))
#     ions_potentials.append(potential(res_ab, res_pos[i], ions_charges, res_pos, current_ion=i))
# print(ions_potentials)

# plt.scatter(x=res_pos/mu, y=ions_potentials, color='r', label=f'$40Ca^+$ ions')

# x = np.linspace(res_pos[0]-1e-6, res_pos[-1]+1e-6)
# curve = res_ab[0] * x**2 + res_ab[1] * x**4
# curve_mean = np.mean(curve)
# # plt.plot(x / mu, curve, label='Trap potential')
# plt.ylabel('V (V)')
# plt.xlabel('x (µm)')
# y = np.ones(len(res_pos)) * curve_mean
# # plt.scatter(x=res_pos/mu, y=y, color='r', label=f'$40Ca^+$ ions')
# plt.legend()
# plt.show()
