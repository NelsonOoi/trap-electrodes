from trapsim import *

valid_positions = [str(i) for i in range(3, 10)]


'''
start: starting electrode index.
end: ending electrode index.
valid_positions: valid positions (labeled by electrode index)
                 for the ion to start and end at.
                 Edge electrodes are not valid positions as they
                 cannot support stable axial confinement.
electrode_width: 120µm by default.
electrode_gap: 5µm by default.

Assumes that the ion is currently centered on the starting electrode position.

'''
s = trap()
'''Configure electrode parameters.'''
l = 1e-6 # µm length scale
u = 20. # V rf peak voltage
m = 40*ct.atomic_mass # 40Ca+ ion mass
q = 1*ct.elementary_charge # ion charge
f = 35e6 #trap rf
o = 2*np.pi*f # rf frequency in rad/s
s['r'].rf = u*np.sqrt(q/m)/(2*l*o)

def axial_shuttle(
        s=s,
        start_pos=[0., 0., 51.7],
        end_pos = [-250., 0., 51.7],
        start='6',
        end='5',
        path_resolution_bound=5,
        electrode_width=120,
        electrode_gap=5,
        valid_positions=valid_positions,
        start_dc_set_file='Apr16_approx_Vs_axial.csv'
    ):
    assert start in valid_positions, 'Ensure that start position is valid.'
    assert end in valid_positions, 'Ensure that end position is valid.'
    assert start != end, 'Ensure that the start and end are not the same.'
    total_displacement = (list(valid_positions).index(end) - list(valid_positions).index(start)) * (electrode_gap + electrode_width)
    keyframes = int(np.ceil(np.abs(total_displacement / path_resolution_bound))) + 1
    print('Total distance: ', total_displacement)
    print('Keyframes: ', keyframes)
    is_going_left = (total_displacement < 0)

    start_dc_set = read_electrode_voltages([start_dc_set_file])[0]
    with (s.with_voltages(start_dc_set)):
        start_pos = s.minimum([0., 0., 51.7])
        start_pos[np.abs(start_pos)<1e-6] = 0

    # obtained from run_trap.py
    # dx2 = 1.6366598620507583e-05
    dx2 = 2e-6
    target_axial_coeffs = [0., dx2, 0., -dx2/2, 0., -dx2/2]
    steps = np.linspace(0, total_displacement, keyframes)
    # print('Steps: ', steps)

    # '''
    # Initial groupings.
    # '''
    # cur_pos_groups = [['1', '11']]
    # cur_pos = int(start)
    # for i in range(-1, 2):
    #     side_pos = cur_pos + i
    #     cur_pos_groups.append([f'{side_pos}', f'{side_pos + 10}'])
    # print('Initial grouping:', cur_pos_groups)

    '''
    cur_pos: exact position for ion to be at.
    cur_electrode_pos: which is the 'center electrode' to control.
    '''
    cur_electrode_pos = int(start)
    print('cur_electrode_pos: ', cur_electrode_pos)
    waveform = []
    electrodes_passed = 0
    perc = 0
    # for step in steps:
    for i in range(0, len(steps)):
        step = steps[i]
        cur_pos = np.array(start_pos) + np.array([step, 0, 0])
        coeffs = get_electrode_coeffs(s=s, ion_pos=cur_pos, save_file=False)
        print('Current position:', cur_pos)
        coeff_indices = np.arange(5)

        '''
        INSERT code to figure out when to swap groupings.
        As the ion moves further away, the further electrodes
        begin to require larger coefficients to make a contribution
        so there has to be a truncation point.

        1. First figure out when to add groups.
        '''
        perc = (np.abs(step) % (electrode_width + electrode_gap)) / (electrode_width + electrode_gap)
        print('perc:', perc)
        print('step:', step)
        print(step % (electrode_width + electrode_gap))
        if (perc > 0.5 and has_passed_midpoint == False):
            '''
            run this once, only after the midpoint
            bewteen each two electrodes has been passed
            '''
            has_passed_midpoint = True
            cur_electrode_pos += np.sign(total_displacement)
            print('Cur pos:', cur_electrode_pos)
        elif (perc <= 0.5):
            has_passed_midpoint = False
        
        '''
        Re-compute groupings based on current position.
        '''
        cur_pos_groups = [['1', '11']]
        for i in range(-1, 2):
            side_electrode_pos = cur_electrode_pos + i
            cur_pos_groups.append([f'{side_electrode_pos}', f'{side_electrode_pos + 10}'])
        print('Initial grouping:', cur_pos_groups)

        electrode_v, group_v = solve_voltages(el_names=s.names, fitted_coeffs=coeffs,
                        target_coeffs=target_axial_coeffs,
                        groups=cur_pos_groups, save_file=False,
                        filename='', coeff_indices=coeff_indices)
        print('Electrode voltages:', electrode_v)
        waveform.append(electrode_v)
    return waveform




print(type(valid_positions))
waveform = axial_shuttle(start='6', end='4', path_resolution_bound=1.)

# def plot_waveform(s, waveform, length=500., res=10001):
#     plt.ion()

#     fig, ax = plt.subplots(1, 1, figsize=(10,7))
#     x0 = s.minimum([0., 0., 51.7])
#     shift = {'z': x0[2]}
#     tmp_len = 20.
#     x3d, x = single_axis('x', (-length/2, length/2), res, shift=shift)
#     for electrode_v in waveform:
#         vx, line1 = 0, 0
#         with s.with_voltages(electrode_v):
#             line1, = ax.plot(x, vx[0], label='Electrode contribution')
#             vx = s.electrical_potential(x3d).T
#     line1.set_ydata(x[0]) 
#     fig.canvas.draw() 
#     fig.canvas.flush_events() 
#     ax.set_title('x-axis')
#     ax.grid()
#     plt.show()

def interactive_plot_waveform(s, waveform, length=1000., res=10001):
    plt.ion()
    fig, ax = plt.subplots(1, 1, figsize=(10,7))
    ax.set_title('x-axis')
    ax.grid()
    ax.set_ylim(-0.08, 0.08)
    x0 = s.minimum([0., 0., 51.7])
    shift = {'z': x0[2]}
    tmp_len = 20.
    x3d, x = single_axis('x', (-length/2, length/2), res, shift=shift)

    electrode_v = waveform[0]
    line1 = 0
    vx = 0
    with s.with_voltages(electrode_v):
        vx = s.electrical_potential(x3d).T
        line1, = ax.plot(x, vx[0], label='Electrode contribution')

    for electrode_v in waveform[1:]:
        with s.with_voltages(electrode_v):
            vx = s.electrical_potential(x3d).T
            line1.set_ydata(vx[0])
        fig.canvas.draw()
        time.sleep(0.03)
        fig.canvas.flush_events() 
    plt.ioff()
interactive_plot_waveform(s=s, waveform=waveform)