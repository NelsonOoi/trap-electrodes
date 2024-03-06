from trapsim import *

electrode_mapping = {
    '1': 10,
    '2': 9,
    '3': 5,
    '4': 8,
    '5': 1,
    '6': 4,
    '7': 3,
    '8': 7,
    '9': 'gnd',
    '10': 6,
    '11': 12,
    '12': 2,
    '13': 11,
    '14': 'gnd',
    '15': 'gnd',
    '16': 'r',
    '17': 'gnd',
    '18': 15,
    '19': 20,
    '20': 14,
    '21': 17,
    '22': 'gnd',
    '23': 16,
    '24': 18,
    '25': 13,
    '26': 19
}

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
print(electrode_ordering)

''' Load trap from GDS.'''
s, electrodes, electrodes_dict = load_trap(
        filename='single_chip.gds',
        ito_layer=12, electrode_mapping=electrode_mapping,
        electrode_ordering=electrode_ordering, plot=False,
        xlim=(-5000,5000), ylim=(-3000,3000), trap_center=trap_center, buildup=False)

'''Configure electrode parameters.'''
l = 1e-6 # µm length scale
u = 10. # V rf peak voltage
m = 40*ct.atomic_mass # 40Ca+ ion mass
q = 1*ct.elementary_charge # ion charge
f = 30e6 #trap rf
o = 2*np.pi*f # rf frequency in rad/s
s['r'].rf = u*np.sqrt(q/m)/(2*l*o)

'''Analyze RF minimum and initial curvatures.'''
# searches for potential minimum point given the starting point
x0 = s.minimum((0, 0, 1.), axis=(1, 2))
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
r = get_axes_potentials(s=s, shift={'z': ion_height})

'''Solve for curvatures of unit voltage on electrodes.'''
curv_filename='1_gds_quetzal.csv'
fc, residuals = get_electrode_curvatures(s, *r,
                            ion_height=ion_height,
                            filename=curv_filename, plot=False)


fc = load_coeffs(filename=curv_filename)
target_axial_curvature = np.array([0., 4e-6, 0., -2e-6, 0.,-2e-6])
# 5 electrode-a-side grouping
groups = [['1', '11'], ['4', '14'], ['6', '16'], ['8', '18'], ['5', '7'], ['15', '17']]
# 7 electrode-a-side grouping
groups = [['1', '11'], ['4', '14', '8', '18'], ['6', '16'], ['3', '13', '9', '19'], ['5', '7'], ['15', '17']]
# groups = [['1', '11'], ['4', '8'], ['14', '18'], ['5', '6', '7', '15', '16', '17'], ['3', '13'], ['9', '19']]
# groups = [['1'],  ['11'], ['4', '8'], ['14', '18'], ['5', '6', '7'], ['15', '16', '17']] # gives unphysical 10^11 results.

Vs_axial_filename = 'Vs_axial_5.csv'
electrode_v, group_v = solve_voltages(el_names=s.names,
                    fitted_coeffs=fc, target_curvature=target_axial_curvature,
                    groups=groups, filename=Vs_axial_filename)
print(electrode_v)

plot_length = 100
plot_fitted_curvature(s=s, electrode_voltages=electrode_v,
                    target_curvature=target_axial_curvature/2,
                    ion_height=ion_height, length=plot_length,
                    shift={'z': ion_height})

# solve_freqs(s=s, f_rad=3e6, f_split=0.2e6, f_axial=1e6, f_traprf=30e6,
#             m=40.*ct.atomic_mass, q=1.*ct.elementary_charge,
#             l=1e-6, u_dc_axial_ref=target_axial_curvature,
#             dc_axial_set_file=Vs_axial_filename,
#             dc_tilt_set_file=Vs_tilt_filename,
#             do_plot_potential=True)

print('\n\n\n')


'''
Solve for tilt modes.
First fit along tilted axes. 
'''
r_tilt = get_axes_potentials(s=s, tilt_theta=np.pi/4, shift={'z': ion_height})

'''Solve for tilt mode curvatures of unit voltage on electrodes.'''
tilt_curv_filename='1_tilt_gds_quetzal.csv'
tilt_fc, residuals = get_electrode_curvatures(s, *r_tilt,
                            ion_height=ion_height,
                            filename=tilt_curv_filename, plot=False)
tilt_fc = load_coeffs(filename=tilt_curv_filename)
target_tilt_curvature = np.array([0., 0., 0., -2e-6, 0., +2e-6])
tilt_groups = [['1'], ['11'], ['4', '8'], ['5', '6', '7'], ['14', '18'], ['15', '16', '17']]

Vs_tilt_filename = 'Vs_tilt.csv'
tilt_electrode_v, tilt_group_v = solve_voltages(el_names=s.names,
                    fitted_coeffs=tilt_fc, target_curvature=target_tilt_curvature,
                    groups=tilt_groups, filename=Vs_tilt_filename)

print(tilt_group_v)

# say that radial frequencies need to be split by x Hz
# assign half of the splitting: x/2 to each
# u_axial = (o_axial + 2*np.pi*(x/2))**2 * m * l**2 / q
# take the difference between the two u_axial components along y' and z'
# compensate using the tilt voltage

solve_freqs(s=s, f_rad=3.5e6, f_split=1e6, f_axial=1e6, f_traprf=30e6,
            m=40.*ct.atomic_mass, q=1.*ct.elementary_charge,
            l=1e-6, u_dc_axial_ref=target_axial_curvature,
            u_dc_tilt_ref=target_tilt_curvature,
            dc_axial_set_file=Vs_axial_filename,
            dc_tilt_set_file=Vs_tilt_filename,
            do_plot_potential=True)
