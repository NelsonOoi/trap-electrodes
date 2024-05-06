from trapsim import *

'''
Algorithm:

1. Select:
    - interpolation (sin, bezier, param).
    - shuttling time (t_shuttle) between start and finish points.
    - final axial phonon number (n_bar) after cooling runs.
    - axial frequency / well curvature during shuttling. usually keep this constant.
2. Calculate amplitude of initial axial oscillation using nbar * hbar * omega_ax(t) = 1/2 * m_Ca+ * omega_ax(t)^2 * A^2
    [x] done.
3. Solve for ion motion over the complete time interval.
   At each time t:
    - Compute position of the well (x_well) using the interpolation functions.
    - Compute axial voltage set for creating a well at x_well.
        - This can be challenging as it will mean finding a way to calculate the voltage set without working from the previous position.
    - Simulate the axial voltage set to find d/dx(Phi) at the ion's position (x_ion).
4. At the final position, FFT the ion's motion and plot it. Deduce the final axial phonon number from this.
5. ???
6. Profit!

Soon:
- Expand to simulations of multi-ion chains.
- Profile optimizer routine to minimize the final oscillation amplitude.


'''

# def test(f_in):
#     t = np.linspace(0, 1, 101)
#     return f_in(t)

# # print(test(sine_interp))

files = ['Apr22_5elec_RF35.0MHz-trapx1.0MHz.csv', 'Vs_axial_karan.csv', 'Vs_axial125_karan_recentered.csv', 'Apr22_gds_Vs_axial.csv', 'Apr22_7dof_gds_Vs_axial.csv']
overall = read_electrode_voltages(files=files)
voltage_set = overall[0]

def get_current_electrode_groups(
                    start_well_pos, start_electrode_pos,
                    cur_well_pos,
                    electrode_width = 120.,
                    electrode_gap = 5.,
                    el_extent = [-2, 2]):
    cur_displacement = np.array(cur_well_pos) - np.array(start_well_pos)
    electrodes_displacement = np.round(cur_displacement / (electrode_width + electrode_gap))
    cur_electrode_pos = start_electrode_pos + int(electrodes_displacement[0])

    perc = (np.abs(cur_displacement[0]) % (electrode_width + electrode_gap)) / (electrode_width + electrode_gap)
    '''
    Re-compute groupings based on current well position.
    '''
    cur_electrode_groups = [['1', '11']]
    for i in range(el_extent[0], el_extent[1]+1):
        side_electrode_pos = int(cur_electrode_pos + i)
        assert side_electrode_pos > 1, 'Invalid shuttling destination.'
        cur_electrode_groups.append([f'{side_electrode_pos}', f'{side_electrode_pos + 10}'])
    return cur_electrode_groups, cur_electrode_pos, perc

def simulate_excitation(s, f_interp, n_phonon=[0.5], modes=[1e6],
                        t_shuttling_start = 10e-6,
                        t_shuttling = 10e-6,
                        q = ct.elementary_charge,
                        ion_mass = (40. * ct.atomic_mass),
                        ion_height = 51.7,
                        start_electrode_pos=6,
                        end_electrode_pos=5,
                        init_dc_set_filename='Apr16_approx_Vs_axial.csv',
                        electrode_width=120.,
                        electrode_gap=5.
                        ):
    '''
    Calculates ion excitation after transport.
    - s:                    Trap system with defined electrodes.
    - f_interp:             Interpolation function from 0 to 1 over time.
    - phonon_n:             Array of phonon numbers along principal axes: (axial, radial1, radial2)
    - modes:                Oscillation frequencies along principal axes: (axial, radial1, radial2)
    - displacement:         Total ion shuttling displacement.
    '''
    init_dc_set = read_electrode_voltages([init_dc_set_filename])[0]
    # Calculate amplitude of initial excitation.
    omega = 2 * np.pi * np.array(modes)
    nhw = np.array(n_phonon) * ct.Planck * omega[0]
    '''
    nhw = 0.5 * m * w^2 * A^2
    A is the initial amplitude of excitation before shuttling.
    '''
    A0 = np.sqrt(2 * nhw / ion_mass / (omega[0]**2))

    # NOTE: Don't need to create the time-range array, because
    # the IVP solver will adjust to find specific time points
    # that may not match up with those found in the linspace array.

    start_well_pos = [(start_electrode_pos - 6) * (electrode_width + electrode_gap),
                      0, ion_height]

    cur_electrode_pos = start_electrode_pos
    total_displacement = (end_electrode_pos - cur_electrode_pos) * (electrode_width + electrode_gap)
    def axial_efield(x, t, t_shuttling_start=t_shuttling_start, t_shuttling=t_shuttling):
        field_x = []
        if (t < t_shuttling_start):
            with s.with_voltages(dcs=init_dc_set):
                field_x = - s.potential([x * m_to_micron, 0, ion_height], 1).flatten()[0] / micron_to_m
            return field_x, start_well_pos
        elif (t >= t_shuttling_start):
            # NOTE solve the shuttling waveform here.
            scaled_t = (t - t_shuttling_start) / t_shuttling
            cur_well_pos = [f_interp(scaled_t) * total_displacement, 0., ion_height]
            coeffs = get_electrode_coeffs(s=s, ion_pos=cur_well_pos, save_file=False)
            coeff_indices = np.arange(5)
            cur_electrode_groups, cur_electrode_pos, perc = get_current_electrode_groups(cur_well_pos=cur_well_pos,
                                                                start_well_pos=start_well_pos,
                                                                start_electrode_pos=start_electrode_pos,
                                                                electrode_width=electrode_width,
                                                                electrode_gap=electrode_gap)
            # print('cur_electrode_groups', cur_electrode_groups)
            electrode_v, group_v = solve_voltages(el_names=s.names, fitted_coeffs=coeffs,
                        target_coeffs=target_axial_coeffs,
                        groups=cur_electrode_groups, save_file=False,
                        filename='', coeff_indices=coeff_indices)
            with s.with_voltages(dcs=electrode_v):
                field_x = - s.potential([x * m_to_micron, 0, ion_height], 1).flatten()[0] / micron_to_m
            return field_x, cur_well_pos

    def fun(t, u):
        x, v = u
        dxdt = v
        field_x, cur_well_pos = axial_efield(x, t)
        dvdt = (q / ion_mass) * field_x
        # print('potential', potential_x, 'x', x, 'dvdt', dvdt, 'v', v)
        print('time:', t * 1e6, 'ion_pos:', x * m_to_micron, 'well pos:', cur_well_pos[0])
        return [dxdt, dvdt]

    u_0 = [A0[0], 0]
    N = 3000
    # t_spacing = 1e-7
    # tmax = N * t_spacing
    t_max = 30e-6
    t_maxplot = 25e-6
    t_pts = np.linspace(0, t_maxplot, N)

    t_span = (0., t_max)
    relerr = 1.e-7
    abserr = 1.e-7
    # NOTE: Start shuttling after t=10 µs. Stay at end position for >10 µs.
    result = solve_ivp(fun, t_span=t_span, y0=u_0, t_eval=t_pts,
                    rtol=relerr, atol=abserr)
    plt.figure()
    plt.plot(t_pts * 1e6, result.y[0,:] * m_to_micron, label = "Numerical solution")
    # plt.plot([math.sin(t) for t in t_pts], "o", label="Analytical solution")
    plt.xlabel("t (µs)")
    plt.ylabel("x (µm)")
    # plt.xlim(0,40)
    # plt.legend()
    plt.grid()
    plt.title('Ion oscillation at axial frequency.')
    plt.show()
    return A0

s, electrodes, electrodes_dict = load_trap()
'''Configure electrode parameters.'''
l = 1e-6 # µm length scale
u = 26.2 # V rf peak voltage
m = 40*ct.atomic_mass # 40Ca+ ion mass
q = 1*ct.elementary_charge # ion charge
f = 35e6 #trap rf
o = 2*np.pi*f # rf frequency in rad/s
s['r'].rf = u*np.sqrt(q/m)/(2*l*o)
x_guess = [0., 0., 50.]
x0 = s.minimum(x_guess)
ion_height = x0[2]

print(simulate_excitation(s=s, f_interp=sine_interp,
                            n_phonon=[0.5], modes=[1e6],
                            t_shuttling_start = 10e-6,
                            t_shuttling = 10e-6,
                            q = ct.elementary_charge,
                            ion_mass = (40. * ct.atomic_mass),
                            ion_height = 51.7,
                            start_electrode_pos=6,
                            end_electrode_pos=5,
                            init_dc_set_filename='Apr16_approx_Vs_axial.csv',
                            electrode_width=120.,
                            electrode_gap=5.))

# print(get_current_electrode_groups(start_pos=[0., 0., 51.7], start_electrode_pos=6,
#                 cur_pos=[187., 0., 51.7], end_electrode_pos=8))