from trapsim import *

today = date.today()
electrode_ordering = [str(i) for i in range(1, 21)]
list2 = ['r', 'gnd']
electrode_ordering = electrode_ordering + list2
trap_center = [(4920+5040)/2, (5502.5 + 5297.5)/2]
trap_center = (0, 0)

is_approx_trap = False
Vs_axial_filename = f'{today.strftime("%b%d")}_gds_Vs_axial.csv'
Vs_tilt_filename = f'{today.strftime("%b%d")}_gds_Vs_tilt.csv'
coeff_filename=f'{today.strftime("%b%d")}_gds_quetzal.csv'
tilt_coeff_filename=f'{today.strftime("%b%d")}_gds_tilt_quetzal.csv'
# tilt_coeff_filename = '1_tilt_gds_quetzal.csv'

s = None
if(is_approx_trap):
    s = approx_trap()
    Vs_axial_filename = f'{today.strftime("%b%d")}_approx_Vs_axial.csv'
    Vs_tilt_filename = f'{today.strftime("%b%d")}_approx_Vs_tilt.csv'
    coeff_filename=f'{today.strftime("%b%d")}_approx_quetzal.csv'
    tilt_coeff_filename=f'{today.strftime("%b%d")}_approx_tilt_quetzal.csv'
    ''' Load trap from GDS.'''
else:
    s, electrodes, electrodes_dict = load_trap(
            filename=append_filepath('single_chip_centered.gds'),
            electrode_layer=default_electrode_config.get('electrode_layer'),
            ito_layer=default_electrode_config.get('ito_layer'), electrode_mapping=default_electrode_config,
            electrode_ordering=electrode_ordering, plot=True,
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
ion_pos = s.minimum((0, 0, 1.), axis=(1, 2))
ion_pos[np.abs(ion_pos)<1e-6] = 0
print('Null point found at: ', ion_pos, 'µm')
print()
ion_height = ion_pos[2]

'''
Obtain potential derivatives due to unit potentials on each electrode.
Method 1. fit second-order polynomial.
'''
# s['r'].rf = 0.
# r = get_axes_unitv_potentials(s=s, shift={'z': ion_height})
# s['r'].rf = u*np.sqrt(q/m)/(2*l*o)
# fc, residuals = get_electrode_coeffs_fit(s, *r,
#                             ion_height=ion_height,
#                             filename=coeff_filename, plot=False)
'''
Method 2. extract using electrode derivative calculator.
'''
derivs = get_electrode_coeffs(s=s, ion_pos=ion_pos, filename=append_filepath(coeff_filename))

'''
Load saved coefficients from file and prepare to solve voltages.
'''
fc = load_coeffs(filename=append_filepath(coeff_filename))
target_axial_coeffs = np.array([0., 2e-6, 0., -1e-6, 0.,-1e-6])

'''
coeff_indices:
indices of which elements from target curvature to use during fit.
'''
coeff_indices = np.arange(5)
# ignore the 2nd order z constraint
# it is redundant due to Gauss's Law imposing a constraint.

# 3 electrode-a-side grouping
# groups = [['1'], ['11'], ['6'], ['16'], ['5', '15'], ['7', '17']] # gives 10^11 results.
# groups = [['1'], ['11'], ['6'], ['16'], ['5', '7'], ['15', '17']] # gives 10^11 results.
# groups = [['1', '11'], ['6', '16'], ['5'], ['15'], ['7'], ['17']] # gives 10^11 results.
# groups = [['6'], ['16'], ['5'], ['7'], ['15'], ['17']] # gives 10^12 results.
# 5 electrode-a-side grouping
# groups = [['1', '11'], ['4', '14'], ['6', '16'], ['8', '18'], ['5', '7'], ['15', '17']] # this one gives physical results within 10V constraint.
# groups = [['1'], ['11'], ['5', '6', '7'], ['15', '16', '17'], ['8', '18'], ['4', '14']]
# 7 electrode-a-side grouping
# groups = [['1', '11'], ['4', '14', '8', '18'], ['6', '16'], ['3', '13', '9', '19'], ['5', '7'], ['15', '17']]
# groups = [['1', '11'], ['4', '8'], ['14', '18'], ['5', '6', '7', '15', '16', '17'], ['3', '13'], ['9', '19']]
# groups = [['1'],  ['11'], ['4', '8'], ['14', '18'], ['5', '6', '7'], ['15', '16', '17']] # gives unphysical 10^11 results.

# 7 DOF case
# groups = [['1', '11'], ['6', '16', '5','15', '7', '17'], ['4', '14', '8', '18', '3', '13', '9', '19']]
# Vs_axial_filename = f'{today.strftime("%b%d")}_7elec_gds_Vs_axial.csv'

# 5 DOF case.
groups = [['1', '11'], ['6', '16', '5','15', '7', '17'], ['4', '14', '8', '18']]
Vs_axial_filename = f'{today.strftime("%b%d")}_5elec_gds_Vs_axial.csv'

# 3 DOF case.
# target_axial_coeffs = np.array([2e-6, -1e-6, -1e-6])
# coeff_indices = [1, 3, 5]
# groups = [['1', '11'], ['6', '16'], ['5', '7', '15', '17']]
# groups = [['1', '11'], ['6', '16'], ['5','15'], ['7', '17']]

print('Using electrode grouping:', groups)
print()

electrode_v, group_v = solve_voltages(el_names=s.names,
                    fitted_coeffs=fc, target_coeffs=target_axial_coeffs,
                    groups=groups, filename=Vs_axial_filename,
                    coeff_indices=coeff_indices, save_file=True)
# with s.with_voltages(electrode_v):
#     print(s.electrical_potential(ion_pos, 'dc', derivative=2, expand=True))
#     print(s.electrical_potential(ion_pos, 'dc', derivative=2, expand=False))

plot_length = 50.
plot_fitted_coeffs(s=s, electrode_voltages=electrode_v,
                    target_coeffs=target_axial_coeffs,
                    ion_height=ion_height, length=plot_length,
                    shift={'z': ion_height},
                    # target_coeff_indices=coeff_indices,
                    plot_target=True)
# plt.show()
# solve_freqs(s=s, f_rad=3e6, f_split=0.2e6, f_axial=1e6, f_traprf=30e6,
#             m=40.*ct.atomic_mass, q=1.*ct.elementary_charge,
#             l=1e-6, u_dc_axial_ref=target_axial_curvature,
#             dc_axial_set_file=Vs_axial_filename,
#             # dc_tilt_set_file=Vs_tilt_filename,
#             do_plot_potential=True)

print('\n\n\n')


'''
Solve for tilt modes.
First fit along tilted axes. 
'''
r_tilt = get_axes_unitv_potentials(s=s, tilt_theta=np.pi/4, shift={'z': ion_height})

'''Solve for tilt curvatures due to unit voltage on each electrode.'''
tilt_fc, residuals = get_electrode_coeffs_fit(s, *r_tilt,
                            ion_height=ion_height,
                            filename=tilt_coeff_filename, plot=False)
tilt_fc = load_coeffs(filename=tilt_coeff_filename)
target_tilt_coeffs = np.array([0., 0., 0., -1e-5, 0., +1e-5])
tilt_groups = [['1'], ['11'], ['5', '7'], ['6'], ['15', '17'], ['16']]
tilt_coeff_indices = np.arange(6)
# tilt_coeff_indices = [0, 1, 2, 3, 4, 5]
# tilt_coeff_indices = [1, 2, 3, 5]
tilt_electrode_v, tilt_group_v = solve_voltages(el_names=s.names,
                    fitted_coeffs=tilt_fc, target_coeffs=target_tilt_coeffs,
                    coeff_indices=tilt_coeff_indices,
                    groups=tilt_groups, filename=Vs_tilt_filename,
                    # exact=True,
                    exact=False)

plot_length = 20.
plot_fitted_coeffs(s=s, electrode_voltages=tilt_electrode_v,
                    target_coeffs=target_tilt_coeffs,
                    ion_height=ion_height, length=plot_length,
                    shift={'z': ion_height},
                    # target_coeff_indices=coeff_indices,
                    plot_target=False)
# with s.with_voltages(tilt_electrode_v):
#     print(s.electrical_potential(ion_pos, 'dc', derivative=2, expand=True))
#     print(s.electrical_potential(ion_pos, 'dc', derivative=2, expand=False))
# plt.show()

'''
say that radial frequencies need to be split by x Hz
assign half of the splitting: x/2 to each
u_axial = (o_axial + 2*np.pi*(x/2))**2 * m * l**2 / q
take the difference between the two u_axial components along y' and z'
compensate using the tilt voltage
'''

solve_freqs_old(s=s, f_rad=3.5e6, f_split=1e6, f_axial=1e6, f_traprf=35e6,
            m=40.*ct.atomic_mass, q=1.*ct.elementary_charge,
            l=1e-6, dc_axial_ref_coeffs=target_axial_coeffs,
            dc_tilt_ref_coeffs=target_tilt_coeffs,
            dc_axial_set_file=Vs_axial_filename,
            dc_tilt_set_file=Vs_tilt_filename,
            do_plot_potential=True,
            save_result=True)
