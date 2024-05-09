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
s = approx_trap()
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
        init_dc_set_filename='Apr16_approx_Vs_axial.csv',
        interp_type='bezier'
    ):
    assert start in valid_positions, 'Ensure that start position is valid.'
    assert end in valid_positions, 'Ensure that end position is valid.'
    assert start != end, 'Ensure that the start and end are not the same.'
    total_displacement = (list(valid_positions).index(end) - list(valid_positions).index(start)) * (electrode_gap + electrode_width)
    keyframes = int(np.ceil(np.abs(total_displacement / path_resolution_bound))) + 1
    print('Total distance: ', total_displacement)
    print('Keyframes: ', keyframes)
    is_going_left = (total_displacement < 0)

    init_dc_set = read_electrode_voltages([init_dc_set_filename])[0]
    with (s.with_voltages(init_dc_set)):
        start_pos = s.minimum(start_pos)
        start_pos[np.abs(start_pos)<1e-6] = 0

    # obtained from run_trap.py
    # dx2 = 1.6366598620507583e-05
    dx2 = 2e-6
    target_axial_coeffs = [0., dx2, 0., -dx2/2, 0., -dx2/2]
    steps = []
    match interp_type:
        case 'linear':
            steps = np.linspace(0, total_displacement, keyframes)
        case 'bezier':
            steps = bezier_steps(keyframes, total_displacement)
        case 'param':
            steps = parametric_steps(keyframes, total_displacement)

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
        left = -2
        right = 2
        for i in range(left, right+1):
            side_electrode_pos = cur_electrode_pos + i
            cur_pos_groups.append([f'{side_electrode_pos}', f'{side_electrode_pos + 10}'])
        print('Initial grouping:', cur_pos_groups)

        electrode_v, group_v = solve_voltages(el_names=s.names, fitted_coeffs=coeffs,
                        target_coeffs=target_axial_coeffs,
                        groups=cur_pos_groups, save_file=False,
                        filename='', coeff_indices=coeff_indices)
        print('Electrode voltages:', electrode_v)
        waveform.append(electrode_v)
    return np.array(waveform)

def bezier_steps(keyframes, total_displacement):
    t = np.linspace(0, 1, keyframes)
    return t * t * (3. - 2. * t) * total_displacement

def parametric_steps(keyframes, total_displacement, alpha=2.):
    t = np.linspace(0, 1, keyframes)
    sq = t**2
    return sq / (alpha * (sq - t) + 1.) * total_displacement

print(type(valid_positions))
bezier_waveform = axial_shuttle(start='6', end='5', path_resolution_bound=1., interp_type='bezier')
# linear_waveform = axial_shuttle(start='6', end='4', path_resolution_bound=1., interp_type='linear')
# param_waveform = axial_shuttle(start='6', end='4', path_resolution_bound=1., interp_type='param')

# waveforms = [bezier_waveform, linear_waveform, param_waveform]
def plot_waveform(s, waveform, length=500., res=10001):
    plt.ioff()
    fig, ax = plt.subplots(1, 1, figsize=(10,7))
    x0 = s.minimum([0., 0., 51.7])
    shift = {'z': x0[2]}
    tmp_len = 20.
    x3d, x = single_axis('x', (-length/2, length/2), res, shift=shift)
    for electrode_v in waveform:
        vx, line1 = 0, 0
        with s.with_voltages(electrode_v):
            vx = s.electrical_potential(x3d).T
            line1, = ax.plot(x, vx[0], label='Electrode contribution')
    line1.set_ydata(x[0]) 
    fig.canvas.draw() 
    fig.canvas.flush_events() 
    ax.set_xlabel('x-axis')
    ax.set_ylabel('Potential (V)')
    ax.set_title('Bezier shuttling profile')
    ax.grid()
    plt.show()

plot_waveform(s, waveform=bezier_waveform)

def show_shuttling_electrode_voltages(s, waveform, title):
    plt.ioff()
    fig, ax = plt.subplots(1, 1, figsize=(10,7))
    ax.grid()
    ax.set_title(title)
    ax.set_ylabel('Voltage (V)')
    ax.set_xlabel('Shuttling timestep')
    for i, name in enumerate(s.names):
        ax.plot(waveform[:, i].T, label=f'El. {name}')
    plt.legend()
    plt.show()

show_shuttling_electrode_voltages(s, bezier_waveform, title='Bezier')

def interactive_plot_compare_waveform(s, waveforms, names, length=1000., res=10001):
    plt.ion()
    fig, ax = plt.subplots(1, 1, figsize=(10,7))
    ax.set_title('x-axis')
    ax.grid()
    ax.set_ylim(-0.08, 0.08)
    x0 = s.minimum([0., 0., 51.7])
    shift = {'z': x0[2]}
    x3d, x = single_axis('x', (-length/2, length/2), res, shift=shift)

    '''
    there is probably a neater way to do this. found it.
    '''
    lines = []
    for i in range(len(waveforms)):
        line1 = 0
        with s.with_voltages(waveforms[i][0]):
            vx = s.electrical_potential(x3d).T
            line1, = ax.plot(x, vx[0], label=f'{names[i]}')
        lines.append(line1)
    plt.legend()
    
    for i in range(1, len(waveforms[0])):
        for j in range(len(waveforms)):
            with s.with_voltages(waveforms[j][i]):
                vx = s.electrical_potential(x3d).T
                lines[j].set_ydata(vx[0])
        fig.canvas.draw()
        # time.sleep(0.03)
        fig.canvas.flush_events()
    plt.ioff()

# interactive_plot_compare_waveform(s=s, waveforms=waveforms, names=['bezier', 'linear', 'parametric'])

def interactive_plot_waveform(s, waveform, length=1000., res=10001):
    '''
    Observe evolution of potential waveforms.
    '''
    plt.ion()
    fig, ax = plt.subplots(1, 1, figsize=(10,7))
    ax.set_title('x-axis')
    ax.grid()
    ax.set_ylim(-0.08, 0.08)
    x0 = s.minimum([0., 0., 51.7])
    shift = {'z': x0[2]}
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
# interactive_plot_waveform(s=s, waveform=waveform)

def demo_interpolations(keyframes=101):
    '''
    Show different interpolation methods.
    '''
    demo_bezier = bezier_steps(keyframes, 1)
    demo_param = parametric_steps(keyframes, 1)
    demo_linear = np.linspace(0, 1, keyframes)
    x = np.arange(keyframes)

    fig, ax = plt.subplots(1, 1, figsize=(7,7))
    # plt.axis('square')
    ax.plot(x, demo_bezier, label='bezier')
    ax.plot(x, demo_linear, label='linear')
    ax.plot(x, demo_param, label='parametric')
    ax.set_ylabel('Progress to completion (%)')
    ax.set_xlabel('Timesteps')
    ax.grid()
    plt.legend()
    plt.ioff()
    plt.show()

# demo_interpolations()