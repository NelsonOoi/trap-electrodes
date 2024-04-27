from trapsim import *

''' Plot trap using GDS '''
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

s = None
approx_trap = False
if (approx_trap):
    s = trap()
else:
    ''' Load trap from GDS.'''
    s, electrodes, electrodes_dict = load_trap(
            filename='single_chip.gds',
            ito_layer=12, electrode_mapping=electrode_mapping,
            electrode_ordering=electrode_ordering, plot=False,
            xlim=(-5000,5000), ylim=(-3000,3000), trap_center=trap_center, buildup=False)
    
length = 6000.
res = 10001
shift = {'z': 51.7}
x3d, x = single_axis('x', bounds=(-length/2,length/2), res=res, shift=shift)


# files = ['Vs_axial_karan.csv', 'Vs_axial125_karan_recentered.csv']
# names = ['Vs_axial', 'Vs_axial125 recentered on electrode 6']

files = ['Vs_axial_karan.csv', 'Vs_axial125_karan_recentered.csv', 'Apr22_gds_Vs_axial.csv', 'Apr22_7dof_gds_Vs_axial.csv']
names = ['Vs_axial', 'Vs_axial125 recentered on electrode 6', 'Apr 22 new 5 electrodes a side', 'Apr 22 new 7 electrodes a side']

# files = ['Vs_axial_karan.csv', 'Apr22_gds_Vs_axial.csv']
# names = ['Vs_axial', 'Apr 22 new 5 electrodes a side']
colors = ['orange', 'green', 'black', 'blue']
# colors = ['orange', 'black']
overall = read_electrode_voltages(files=files)

plt.figure(figsize=(10,5))
for i in range(len(overall)):
    voltage_set = overall[i]
    print(files[i])
    fit_order = 2
    print(voltage_set)
    with s.with_voltages(dcs=voltage_set):
        p_x = s.electrical_potential(x3d, typ='dc', derivative=0)
        plt.plot(x, p_x, color=colors[i], label=names[i])
plt.ylabel('Potential (V)')
plt.xlabel('x-axis (µm)')
plt.legend()
plt.grid()
plt.show()