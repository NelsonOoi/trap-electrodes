from trapsim import *
prefix = 'Electrode'
suffix = 'line_xion0[um]_yion0 [um]0_zion51.8 [um].csv'

n_electrodes = 10
el_nums = np.arange(1, n_electrodes + 1)
axis_names = ['x', 'y', 'z']

s, electrodes, electrodes_dict = load_trap()
'''Configure electrode parameters.'''
l = 1e-6 # µm length scale
u = 20. # V rf peak voltage
m = 40*ct.atomic_mass # 40Ca+ ion mass
q = 1*ct.elementary_charge # ion charge
f = 35e6 #trap rf
o = 2*np.pi*f # rf frequency in rad/s
s['r'].rf = u*np.sqrt(q/m)/(2*l*o)

ion_height = 51.8
shift = {'z': ion_height}
res = 10001

for el in el_nums:
    fig, ax = plt.subplots(1, 3, figsize=(20,6))
    fig.suptitle(f'Electrode {el}.')
    ax = ax.flatten()
    for i, axis in enumerate(axis_names): 
        filename = prefix + str(el) + '_' + axis + suffix
        filename = append_filepath(filename='') + 'COMSOL/' + filename
        # print(filename)
        if (axis == 'x'):
            el_pot = pd.read_csv(filename, skiprows=8, header=None)
        else:
            el_pot = pd.read_csv(filename, usecols=[axis, 'V2 (V)'], skiprows=8)
        # print(el_pot)
        el_pot = el_pot.sort_values(by=el_pot.columns[0])
        
        ax[i].plot(el_pot.iloc[:, 0], el_pot.iloc[:, 1], label=f'COMSOL')
        ax[i].set_xlabel(f'{axis} (µm)')
        ax[i].grid()

        '''Prepare to plot Python simulation.'''
        a3d, a = single_axis(axis, bounds=(el_pot[el_pot.columns[0]].min(),
                                           el_pot[el_pot.columns[0]].max()), res=res, shift=shift)
        print(a3d)
        unit_voltage = np.zeros(len(s.names))
        unit_voltage[el-1] = 1.
        with s.with_voltages(unit_voltage):
            p_a = s.electrical_potential(a3d, 'dc', 0)
            print(p_a.flatten())
            ax[i].plot(a, p_a.flatten(), label=f'Python')
            ax[i].legend()

plt.show()
    
    # ax.plot
    

