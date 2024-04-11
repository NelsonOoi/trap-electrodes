from trapsim import *
import matplotlib as mpl

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

# fig, ax = plt.subplots(1, 2, figsize=(15,10))
# print(s["2"].to_points().points) # use this to get the centroid coordinates of an electrode
# s.plot(ax[0])
# # s.plot_voltages(ax[1], u=s.rfs)
# # s['2'].dc = 1.
# s.plot_voltages(ax[1])
# r= 1000
# for axi in ax.flat:
#     axi.set_aspect("equal")
#     axi.set_xlim(-r, r)
#     axi.set_ylim(-r, r)



fig, ax = plt.subplots(1, 2, figsize=(15,5))
# print(s["2"].to_points().points) # use this to get the centroid coordinates of an electrode
s.plot(ax[0])
# s.plot_voltages(ax[1], u=s.rfs)
# s['2'].dc = 1.
s.plot_voltages(ax[1])
r= 2000
for axi in ax.flat:
    axi.set_aspect("equal")
    axi.set_xlim(-r, r)
    axi.set_ylim(-r, r)

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
x0 = s.minimum((0., 0., 1.), axis=(1, 2))
ion_height = x0[2]
ion_height = 51.7
print('x0', x0, 'ion height', ion_height)

print(s.names)

# '''Extract alphas and betas from voltage sets.'''
# files = ['Vs_axial_karan.csv', 'Vs_axial125_karan_recentered.csv', 'Vs_tilt_karan.csv', 'Vs_xcomp.csv', 'Vs_ycomp.csv', 'Vs_zcomp.csv']
# names = ['Vs_axial.', 'Vs_axial125 recentered at electrode 6.', 'Vs_tilt.', 'Vs_xcomp.', 'Vs_ycomp.', 'Vs_zcomp.']
# files=['Vs_axial_karan2.csv']
files=['apr5_gds_Vs_axial.csv']
names=['adskhaksdjha']

overall = read_electrode_voltages(files=files)
print(overall)
for i in range(len(overall)):
    voltage_set = overall[i]
    print(files[i])
    fit_order = 2
    print(voltage_set)
    with s.with_voltages(dcs=voltage_set):

        '''Temporary testing code'''
        fig1, ax1 = plt.subplots(1, 2, figsize=(15,5))
        s.plot(ax1[0])
        s.plot_voltages(ax1[1])
        r= 2000
        for axi in ax1.flat:
            axi.set_aspect("equal")
            axi.set_xlim(-r, r)
            axi.set_ylim(-r, r)
        
        u = s.dcs
        um = np.fabs(u).max() or 1.
        u = (u / um + 1)/2
        cmap = plt.cm.RdBu_r
        print('min u:', np.min(u), 'max u:', np.max(u), 'u', u)
        # norm = mpl.colors.Normalize(vmin=np.min(u), vmax=np.max(u))
        norm = mpl.colors.Normalize(vmin=0, vmax=1)
        # colors = [cmap(ui) for ui in u]
        # print(colors)
        # bar_range = np.round(np.arange(-um, um, 2*um/10), 3)
        # code to make the colorbar
        bar_range = np.round(np.linspace(-um, um, 11), 3)
        colorb = fig1.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=ax1, orientation='vertical', label='Applied voltages (V)')
        colorb.ax.yaxis.set_major_locator(plt.MaxNLocator(len(bar_range)))
        colorb.ax.set_yticklabels(bar_range)
        # fig1.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        #      ax=ax1, orientation='vertical', label='Applied voltages (V)').ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base=0.1))
        '''temp test code end'''

        r1 = s.electrical_potential(x0, typ="dc", derivative=1, expand=True)
        r2 = s.electrical_potential(x0, typ="dc", derivative=2, expand=True)
        print('e-field:', r1)
        print('curvatures:', r2)
        r = get_axes_potentials(s=s, length=3000., shift={'z': ion_height})
        fig, ax = plt.subplots(1, 3, figsize=(10,5))
        ax = ax.flatten()
        axisname = 'xyz'
        for j in range(3):
            axes = r[j]
            pots = r[j+3].flatten()
            print(len(axes), len(pots))
            
            ax[j].plot(axes, pots)
            ax[j].set_xlabel(f'{axisname[j]} (µm)')
        ax[0].set_ylabel('Potential (V)')
        fig.suptitle(f'{names[i]}')
        print('voltage set', s.electrical_potential([0., 100., 1.], derivative=0))
plt.show()
