from trapsim import *

# used to center figure
center_electrodes = ['6', '16']
trap_center = [(4920+5040)/2, (5502.5 + 5297.5)/2]

alternate_ordering = False
electrode_ordering = [str(i) for i in range(1, 21)]
list2 = ['r', 'gnd']

if (alternate_ordering):
    electrode_ordering = [str(i) for i in range(2, 21)]
    electrode_ordering.remove('11')
    list2 = ['1', '11', 'r', 'gnd']

electrode_ordering = electrode_ordering + list2

coeff_filename = ''
approx_trap = False
s = None
if(approx_trap):
    s=trap()
    coeff_filename='apr9_approx_quetzal.csv'
else:
    coeff_filename='apr9_electrode_coeffs_gds_quetzal.csv'
    # coeff_filename = 'apr9_fitted_curvature_gds_quetzal.csv'
    # coeff_filename = 'apr9_test_fullstop_fitted_curvature_gds_quetzal.csv'
    ''' Load trap from GDS.'''
    s, electrodes, electrodes_dict = load_trap(
            filename='single_chip.gds',
            electrode_layer=37,
            ito_layer=12, electrode_mapping=default_electrode_mapping,
            electrode_ordering=electrode_ordering, plot=False,
            xlim=(-5000,5000), ylim=(-3000,3000), trap_center=trap_center, buildup=False)

'''Configure electrode parameters.'''
l = 1e-6 # µm length scale
u = 20. # V rf peak voltage
m = 40*ct.atomic_mass # 40Ca+ ion mass
q = 1*ct.elementary_charge # ion charge
f = 35e6 #trap rf
o = 2*np.pi*f # rf frequency in rad/s
s['r'].rf = u*np.sqrt(q/m)/(2*l*o)

'''Analyze RF minimum and initial curvatures.'''
# searches for potential minimum point given the starting point
x0 = s.minimum((0, 0, 1.), axis=(1, 2))
x0[np.abs(x0)<1e-6] = 0
print('x0', x0)
ion_pos = x0
ion_height = x0[2]

''' Plot trap using electrode and verify middle electrode position'''
# vs=np.zeros(len(s.names))
# vs[4] = 1
# with (s.with_voltages(vs)):
#     fig, ax = plt.subplots(1, 2, figsize=(15,10))
#     s.plot(ax[0])
#     s.plot_voltages(ax[1])
#     r=3000
#     for axi in ax.flat:
#         axi.set_aspect('equal')
#         axi.set_xlim(-r, r)
#         axi.set_ylim(-r, r)
#     plt.show()
# for line in s.analyze_static(x=x0, axis=(1, 2), m=m, q=q, l=l, o=o):
#     print(line)

'''Obtain curvatures due to unit potentials on each electrode.'''
# s['r'].rf = 0.
# r = get_axes_unitv_potentials(s=s, shift={'z': ion_height})
# s['r'].rf = u*np.sqrt(q/m)/(2*l*o)
'''Solve for curvatures of unit voltage on electrodes.'''

derivs = get_electrode_coeffs(s=s, ion_pos=ion_pos, filename=coeff_filename)


# fc, residuals = get_electrode_curvatures(s, *r,
#                             ion_height=ion_height,
#                             filename=coeff_filename, plot=False)


fc = load_coeffs(filename=coeff_filename)


target_axial_curvature = np.array([0., 4e-6, 0., -2e-6, 0.,-2e-6])
coeff_indices = np.arange(6)
# target_axial_curvature = np.array([0., 4e6, 0., -2e6, 0.,-2e6])
# 3 electrode-a-side grouping
# groups = [['1'], ['11'], ['6'], ['16'], ['5', '15'], ['7', '17']] # gives 10^11 results.
# groups = [['1'], ['11'], ['6'], ['16'], ['5', '7'], ['15', '17']] # gives 10^11 results.
groups = [['1', '11'], ['6', '16'], ['5'], ['15'], ['7'], ['17']] # gives 10^11 results.
# groups = [['6'], ['16'], ['5'], ['7'], ['15'], ['17']] # gives 10^12 results.
# 5 electrode-a-side grouping
# groups = [['1', '11'], ['4', '14'], ['6', '16'], ['8', '18'], ['5', '7'], ['15', '17']] # this one gives physical results within 10V constraint.
# groups = [['1'], ['11'], ['5', '6', '7'], ['15', '16', '17'], ['8', '18'], ['4', '14']]
# 7 electrode-a-side grouping
# groups = [['1', '11'], ['4', '14', '8', '18'], ['6', '16'], ['3', '13', '9', '19'], ['5', '7'], ['15', '17']]
# groups = [['1', '11'], ['4', '8'], ['14', '18'], ['5', '6', '7', '15', '16', '17'], ['3', '13'], ['9', '19']]
# groups = [['1'],  ['11'], ['4', '8'], ['14', '18'], ['5', '6', '7'], ['15', '16', '17']] # gives unphysical 10^11 results.

# 5 DOF case.

# 3 DOF case.
# target_axial_curvature = np.array([4e-6, -2e-6, -2e-6])
# groups = [['1', '11'], ['6', '16'], ['5', '7', '15', '17']]
# coeff_indices = [1, 3, 5]

print('Using group:', groups)

Vs_axial_filename = 'apr10_gds_Vs_axial.csv'
electrode_v, group_v = solve_voltages(el_names=s.names,
                    fitted_coeffs=fc, target_curvature=target_axial_curvature,
                    groups=groups, filename=Vs_axial_filename,
                    coeff_indices=coeff_indices
                    )
print(electrode_v)

plot_length = 50.

target_coeff_indices = [1, 3, 5]
if (len(target_axial_curvature) == 3):
    target_coeff_indices = [0, 1, 2]

plot_fitted_curvature(s=s, electrode_voltages=electrode_v,
                    target_curvature=target_axial_curvature,
                    ion_height=ion_height, length=plot_length,
                    shift={'z': ion_height},
                    target_coeff_indices=target_coeff_indices,
                    plot_target=False
                    # plot_target=True
                    )
plt.show()

solve_freqs(s=s, f_rad=3e6, f_split=0.2e6, f_axial=1e6, f_traprf=30e6,
            m=40.*ct.atomic_mass, q=1.*ct.elementary_charge,
            l=1e-6, u_dc_axial_ref=target_axial_curvature,
            dc_axial_set_file=Vs_axial_filename,
            # dc_tilt_set_file=Vs_tilt_filename,
            do_plot_potential=True)

print('\n\n\n')


# '''
# Solve for tilt modes.
# First fit along tilted axes. 
# '''
# r_tilt = get_axes_unitv_potentials(s=s, tilt_theta=np.pi/4, shift={'z': ion_height})

# '''Solve for tilt mode curvatures of unit voltage on electrodes.'''
# tilt_coeff_filename='1_tilt_gds_quetzal.csv'
# tilt_fc, residuals = get_electrode_curvatures(s, *r_tilt,
#                             ion_height=ion_height,
#                             filename=tilt_coeff_filename, plot=False)
# tilt_fc = load_coeffs(filename=tilt_coeff_filename)
# target_tilt_curvature = np.array([0., 0., 0., -2e-6, 0., +2e-6])
# tilt_groups = [['1'], ['11'], ['4', '8'], ['5', '6', '7'], ['14', '18'], ['15', '16', '17']]

# Vs_tilt_filename = 'mar28_Vs_tilt.csv'
# tilt_electrode_v, tilt_group_v = solve_voltages(el_names=s.names,
#                     fitted_coeffs=tilt_fc, target_curvature=target_tilt_curvature,
#                     groups=tilt_groups, filename=Vs_tilt_filename)

# print(tilt_group_v)

# # say that radial frequencies need to be split by x Hz
# # assign half of the splitting: x/2 to each
# # u_axial = (o_axial + 2*np.pi*(x/2))**2 * m * l**2 / q
# # take the difference between the two u_axial components along y' and z'
# # compensate using the tilt voltage

# solve_freqs(s=s, f_rad=3.5e6, f_split=1e6, f_axial=1e6, f_traprf=35e6,
#             m=40.*ct.atomic_mass, q=1.*ct.elementary_charge,
#             l=1e-6, u_dc_axial_ref=target_axial_curvature,
#             u_dc_tilt_ref=target_tilt_curvature,
#             dc_axial_set_file=Vs_axial_filename,
#             dc_tilt_set_file=Vs_tilt_filename,
#             do_plot_potential=True)
