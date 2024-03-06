from trapsim import *
import pandas as pd

'''
Boilerplate for trap loading and configuration.
'''
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

# files = ['Vs_axial_karan.csv', 'Vs_axial125_karan_recentered.csv', 'Vs_tilt_karan.csv', 'Vs_xcomp.csv', 'Vs_ycomp.csv', 'Vs_zcomp.csv']
files = ['Vs_test.csv']
for filename in files:
    filename = append_filepath(filename=filename)
    el_v = pd.read_csv(filename, header=None).to_numpy().flatten()
    el_v = np.append(el_v, np.array([0, 0]))
    print(el_v)
    with s.with_voltages(el_v):
        plot_field(s=s,
                    x_grid_bounds=(-30., 30.),
                    y_grid_bounds=(-10., 10.),
                    z_grid_bounds=(30., 70.),
                #    x_grid_bounds=(-300., 300.),
                #     y_grid_bounds=(5., 500.),
                #     z_grid_bounds=(0.1, 50),
                    grid_res=(101, 101, 101),
                    is_dc_potential=True)

plt.show()

