from trapsim import *

from tkinter import *
from tkinter import filedialog, ttk
from tkinter.messagebox import showinfo
import os
                                                                                             
# Create the root window
window = Tk()
# Set window title
window.title('Electrode Simulation')
# Set window size
window.geometry('900x700')
options = { 'padx': 5, 'pady': 5}

'''
Trap type selecton.
'''
def show_selected_trap_type():
    showinfo(
        title='Result',
        message=selected_trap_type.get()
    )
    print(selected_trap_type.get())
    # if(selected_trap_type.get() == trap_types[1][1]):
        # label_file_explorer.pack(fill='x', padx=5, pady=5)
        # button_explore.pack(fill='x', padx=5, pady=5)

trap_types = (('Approximate Quetzal (Symmetric).', 'Approx'),
         ('GDS', 'GDS'))

# label
trap_label = ttk.Label(text="Select trap type:")
# label.pack(fill='x', padx=5, pady=5)
trap_label.grid(column=0, row=0, **options)

# set default selected trap type.
selected_trap_type = StringVar(value=trap_types[0][1])
# radio buttons
def radio_command():
    if(selected_trap_type.get() == trap_types[1][1]):
        # label_file_explorer.pack(fill='x', padx=5, pady=5)
        # button_explore.pack(fill='x', padx=5, pady=5)
        label_file_selected.grid(column=1, row=2, **options)
        button_explore.grid(column=2, row=2, **options)
    else:
        # label_file_explorer.pack_forget()
        # button_explore.pack_forget()
        label_file_selected.grid_forget()
        button_explore.grid_forget()

for i, trap_type in enumerate(trap_types):
    r = ttk.Radiobutton(
        window,
        text=trap_type[0],
        value=trap_type[1],
        variable=selected_trap_type,
        command=radio_command
    )
    # r.pack(fill='x', padx=5, pady=5)
    r.grid(column=i+1, row=0, **options)

# button
# button = ttk.Button(
#     window,
#     text="Get Selected Trap",
#     command=show_selected_trap_type)
# button.pack(fill='x', padx=5, pady=5)
  
#Set window background color
# window.config(background = "gray")
  
# Function for opening the 
# file explorer window
gds_filename = ''
def browse_GDS_files():
    global gds_filename
    gds_filename = filedialog.askopenfilename(initialdir = os.getcwd(),
                                          title = 'Select a File',
                                          filetypes = (('GDS files',
                                                        '*.gds'),
                                                       ('All files',
                                                        '*.*')))
      
    # Change label contents
    label_file_selected.configure(text="File: " + '.../' + gds_filename.split('/')[-2] + '/' + gds_filename.split('/')[-1])

electrode_config_filename = ''
def browse_electrode_JSON_files():
    global electrode_config_filename
    electrode_config_filename = filedialog.askopenfilename(initialdir = os.getcwd(),
                                          title = 'Select a File',
                                          filetypes = (('JSON files',
                                                        '*.json'),
                                                       ('All files',
                                                        '*.*')))
    label_config_file.configure(text="File: " + '.../' + electrode_config_filename.split('/')[-2] + '/' + electrode_config_filename.split('/')[-1])

# if (selected_trap_type.get() == trap_types[1][1]):
label_config_file = Label(window, 
                            text = 'Default config selected.',
                            height = 1,
                            fg='red')
# label_file_explorer.pack(fill='x', padx=5, pady=5)
    
label_file_selected = Label(window, 
                            text = 'No trap selected.',
                            height = 1,
                            fg='red')

button_explore = ttk.Button(window, 
                        text = 'Select Trap GDS',
                        command = browse_GDS_files)


button_explore_config = ttk.Button(window, 
                        text = 'Select Config JSON',
                        command = browse_electrode_JSON_files)


'''
Input trap parameters.
'''
row = 1

label_config_file.grid(column=1, row=row, **options)
button_explore_config.grid(column=2, row=row, **options)

row += 2
# trap rf - f
rf_label = Label(window, 
                text = 'RF (MHz):',
                height = 1)
rf_text = StringVar(value='35.0')
rf_entry = ttk.Entry(textvariable=rf_text)
rf_label.grid(column=0, row=row, **options)
rf_entry.grid(column=1, row=row, **options)

row += 1
# single ion charge - q
charge_label = Label(window, 
                text = 'Ion charge (+e):',
                height = 1)
charge_text = StringVar(value='1.0')
charge_entry = ttk.Entry(textvariable=charge_text)
charge_label.grid(column=0, row=row, **options)
charge_entry.grid(column=1, row=row, **options)

def artiq_export_check_command(row=row):
    if(is_artiq_export.get() == 1):
        artiq_button_explore.grid(column=2, row=row+1, **options)
        artiq_label_file_selected.grid(column=2, row=row+2, **options)
    else:
        artiq_button_explore.grid_forget()
        artiq_label_file_selected.grid_forget()

artiq_config_filename = ''
def browse_artiq_JSON_files():
    global artiq_config_filename
    artiq_config_filename = filedialog.askopenfilename(initialdir = os.getcwd(),
                                          title = 'Select a File',
                                          filetypes = (('JSON files',
                                                        '*.json'),
                                                       ('All files',
                                                        '*.*')))
    artiq_label_file_selected.configure(text="File: " + '.../' + artiq_config_filename.split('/')[-2] + '/' + artiq_config_filename.split('/')[-1])

is_artiq_export = IntVar(value=1)
artiq_export_check = ttk.Checkbutton(window, text='Export to ARTIQ-readable format?',
                                     variable=is_artiq_export, onvalue=1, offvalue=0,
                                     command=artiq_export_check_command)
artiq_export_check.grid(column=2, row=row, **options)

row += 1
artiq_button_explore = ttk.Button(window, 
                        text = 'ARTIQ DAC map JSON',
                        command = browse_artiq_JSON_files)
artiq_button_explore.grid(column=2, row=row, **options)

# single ion mass - m
mass_label = Label(window, 
                text = 'Ion mass (amu):',
                height = 1)
mass_text = StringVar(value='40.0')
mass_entry = ttk.Entry(textvariable=mass_text)
mass_label.grid(column=0, row=row, **options)
mass_entry.grid(column=1, row=row, **options)


row += 1
# lengthscale - l
lengthscale_label = Label(window, 
                text = 'Trap length scale (µm):',
                height = 1, 
                fg = 'white')
lengthscale_text = StringVar(value='1.0')
lengthscale_entry = ttk.Entry(textvariable=lengthscale_text)
lengthscale_label.grid(column=0, row=row, **options)
lengthscale_entry.grid(column=1, row=row, **options)

artiq_label_file_selected = Label(window, 
                            text = 'Default ARTIQ DAC map selected.',
                            height = 1,
                            fg='red')
artiq_label_file_selected.grid(column=2, row=row, **options)

row += 1
# target axial frequency
x_axial_label = Label(window, 
                text = 'Target axial freq. (MHz):',
                height = 1, 
                fg = 'white')
x_axial_text = StringVar(value='1.0')
x_axial_entry = ttk.Entry(textvariable=x_axial_text)
x_axial_label.grid(column=0, row=row, **options)
x_axial_entry.grid(column=1, row=row, **options)

row += 1
y_radial_label = Label(window, 
                text = 'Target radial y freq. (MHz):',
                height = 1, 
                fg = 'white')
y_radial_text = StringVar(value='3.0')
y_radial_entry = ttk.Entry(textvariable=y_radial_text)
y_radial_label.grid(column=0, row=row, **options)
y_radial_entry.grid(column=1, row=row, **options)

def translated_well_check_command(row=row):
    if(is_translated_export.get() == 1):
        translated_well_entry.grid(column=2, row=row+1, **options)
    else:
        translated_well_entry.grid_forget()

is_translated_export = IntVar(value=1)
translated_well_check = ttk.Checkbutton(window, text='Export off-center axial well? (Enter electrode position.)',
                                     variable=is_translated_export, onvalue=1, offvalue=0,
                                     command=translated_well_check_command)
translated_well_check.grid(column=2, row=row, **options)


row += 1

translated_well_text = StringVar(value='6')
translated_well_entry = ttk.Entry(textvariable=translated_well_text)
translated_well_entry.grid(column=2, row=row, **options)

z_radial_label = Label(window, 
                text = 'Target radial z freq. (MHz):',
                height = 1, 
                fg = 'white')
z_radial_text = StringVar(value='4.0')
z_radial_entry = ttk.Entry(textvariable=z_radial_text)
z_radial_label.grid(column=0, row=row, **options)
z_radial_entry.grid(column=1, row=row, **options)

row += 1
tilt_angle_label = Label(window, 
                text = 'Radial tilt angle (deg):',
                height = 1, 
                fg = 'white')
tilt_angle_text = StringVar(value='45.0')
tilt_angle_entry = ttk.Entry(textvariable=tilt_angle_text)
tilt_angle_label.grid(column=0, row=row, **options)
tilt_angle_entry.grid(column=1, row=row, **options)

row += 1
basis_axial_label = Label(window, 
                text = 'Basis axial set coefficients\n <x1 (V/m), x2 (V/m^2), \ny1 (V/m), y2 (V/m^2),\nz1 (V/m), z2 (V/m^2)>\n comma separated',
                height = 5, 
                fg = 'white')
basis_axial_text = StringVar(value='0, 2e6, 0, -1e6, 0, -1e6')
basis_axial_entry = ttk.Entry(textvariable=basis_axial_text)
basis_axial_label.grid(column=0, row=row, **options)
basis_axial_entry.grid(column=1, row=row, **options)

row += 1
basis_tilt_label = Label(window, 
                text = 'Basis tilt set coefficients\n <x1 (V/m), x2 (V/m^2), \ny1 (V/m), y2 (V/m^2),\nz1 (V/m), z2 (V/m^2)>\n comma separated',
                height = 5, 
                fg = 'white')
basis_tilt_text = StringVar(value='0, 0, 0, -1e7, 0, 1e7')
basis_tilt_entry = ttk.Entry(textvariable=basis_tilt_text)
basis_tilt_label.grid(column=0, row=row, **options)
basis_tilt_entry.grid(column=1, row=row, **options)

# rf_label.pack(fill='x', padx=5, pady=5)
# rf_text.pack(fill='x', padx=5, pady=5)
# button_explore.pack(fill='x', padx=5, pady=5)

# button_print = Button(window, 
#                     text = "Print filename",
#                     command = lambda: print(gds_filename)) 

# button_exit = Button(window, 
#                      text = "Exit",
#                      command = exit) 
# button_exit.pack(fill='x', padx=5, pady=5)

# Grid method is chosen for placing
# the widgets at respective positions 
# in a table like structure by
# specifying rows and columns
# label_file_explorer.grid(column = 1, row = 1)
  
# button_explore.grid(column = 2, row = 1)
  
# button_exit.grid(column = 1,row = 3)

def read_input_coeffs(text):
    target_coeffs = np.float64(text.get().replace(' ', '').split(','))
    for i, coeff in enumerate(target_coeffs):
        target_coeffs[i] *= micron_to_m
        if (i % 2 == 1):
            target_coeffs[i] *= micron_to_m
    return target_coeffs

def read_input_ints(text):
    int_val = int(text.get().replace(' ', ''))
    return int_val

def read_input_floats(text):
    float_val = np.float64(text.get().replace(' ', ''))
    return float_val

def run_simulation():
    today = date.today()
    is_gds = False
    if(selected_trap_type.get() == trap_types[1][1]):
        is_gds = True
    f_traprf = mhz * read_input_floats(rf_text)
    q = ct.elementary_charge * read_input_floats(charge_text)
    m = ct.atomic_mass * read_input_floats(mass_text)
    l = micron_to_m * read_input_floats(lengthscale_text)
    f_axial = mhz * read_input_floats(x_axial_text)
    f_rad_y = mhz * read_input_floats(y_radial_text)
    f_rad_z = mhz * read_input_floats(z_radial_text)
    tilt_angle_radians = np.pi / 180. * read_input_floats(tilt_angle_text)
    target_axial_coeffs = read_input_coeffs(basis_axial_text)
    target_tilt_coeffs = read_input_coeffs(basis_tilt_text)
    translated_well_electrode_position = read_input_ints(translated_well_text)

    '''
    Details specific to simulation.
    '''

    '''
    Create directory for each run.
    '''
    rundir = append_filepath(os.path.join('runs', datetime.now().strftime('%Y-%m-%d-%H-%M-%S')))
    if (not os.path.exists(rundir)):
        os.makedirs(rundir)
    print(rundir)

    u_temp = 26.2 #V
    o_traprf = 2 * np.pi * f_traprf
    print(q, m, l, target_axial_coeffs, target_tilt_coeffs)
    s = approx_trap()

    '''
    Config file must contain:
    - axial set groups.
    - tilt set groups.
    - Only if GDS option selected, must also include:
        - electrode_layer to extract from GDS.
        - ito_layer.
        - electrode polygon name/ordering.
    '''
    electrode_config = default_electrode_config
    if (electrode_config_filename != ''):
        with open(electrode_config_filename) as f:
            electrode_config = json.load(f)
    print(electrode_config)

    translate_well_by = int(translated_well_electrode_position - electrode_config.get('trap_center_electrode'))
    axial_electrodes = electrode_config.get('axial_electrodes')

    artiq_config = default_artiq_config
    if (artiq_config_filename != ''):
        with open(artiq_config_filename) as f:
            artiq_config = json.load(f)
    print(artiq_config)

    Vs_axial_filename = f'{today.strftime("%b%d")}_approx_Vs_axial.csv'
    Vs_tilt_filename = f'{today.strftime("%b%d")}_approx_Vs_tilt.csv'
    Vs_scaledfreq_axial_filename = f'{today.strftime("%b%d")}_approx_Vs_axial_{np.round(f_axial/1e6)}MHz.csv'
    Vs_scaledfreq_tilt_filename = f'{today.strftime("%b%d")}_approx_Vs_tilt_y_{np.round(f_rad_y/1e6)}MHz_z_{np.round(f_rad_z/1e6)}MHz.csv'
    Vs_scaledfreq_overall_filename = f'{today.strftime("%b%d")}_approx_Vs_overall_x_{np.round(f_axial/1e6)}MHz_y_{np.round(f_rad_y/1e6)}MHz_z_{np.round(f_rad_z/1e6)}MHz.csv'
    Vs_artiq_axial_filename = f'{today.strftime("%b%d")}_ARTIQ_approx_Vs_axial_{np.round(f_axial/1e6)}MHz.csv'
    Vs_artiq_tilt_filename = f'{today.strftime("%b%d")}_ARTIQ_approx_Vs_tilt_y_{np.round(f_rad_y/1e6)}MHz_z_{np.round(f_rad_z/1e6)}MHz.csv'

    coeff_filename=f'{today.strftime("%b%d")}_approx_quetzal.csv'
    tilt_coeff_filename=f'{today.strftime("%b%d")}_approx_tilt_quetzal.csv'
    if (is_gds):
        print('is gds', is_gds)
        Vs_axial_filename = f'{today.strftime("%b%d")}_gds_Vs_axial.csv'
        Vs_tilt_filename = f'{today.strftime("%b%d")}_gds_Vs_tilt.csv'
        Vs_scaledfreq_axial_filename = f'{today.strftime("%b%d")}_gds_Vs_axial_{np.round(f_axial/1e6)}MHz.csv'
        Vs_scaledfreq_tilt_filename = f'{today.strftime("%b%d")}_gds_Vs_tilt_y_{np.round(f_rad_y/1e6)}MHz_z_{np.round(f_rad_z/1e6)}MHz.csv'
        Vs_scaledfreq_overall_filename = f'{today.strftime("%b%d")}_gds_Vs_overall_x_{np.round(f_axial/1e6)}MHz_y_{np.round(f_rad_y/1e6)}MHz_z_{np.round(f_rad_z/1e6)}MHz.csv'

        Vs_scaledfreq_artiq_axial_filename = f'{today.strftime("%b%d")}_ARTIQ_gds_Vs_axial_{np.round(f_axial/1e6)}MHz.csv'
        Vs_scaledfreq_artiq_tilt_filename = f'{today.strftime("%b%d")}_ARTIQ_gds_Vs_tilt_y_{np.round(f_rad_y/1e6)}MHz_z_{np.round(f_rad_z/1e6)}MHz.csv'

        coeff_filename=f'{today.strftime("%b%d")}_gds_quetzal.csv'
        tilt_coeff_filename=f'{today.strftime("%b%d")}_gds_tilt_quetzal.csv'
        
        s, electrodes, electrodes_dict = load_trap(filename=gds_filename,
                      electrode_layer=electrode_config.get('electrode_layer'),
                      ito_layer=electrode_config.get('ito_layer'),
                      plot=True, buildup=False, electrode_mapping=electrode_config,
                      trap_center=(0, 0), xlim=(-4000, 4000), ylim=(-3000, 3000))
    s['r'].rf = u_temp*np.sqrt(q/m)/(2*l*o_traprf)
    print(s.names)

    '''Analyze RF minimum and initial curvatures.'''
    # searches for potential minimum point given the starting point
    ion_pos = s.minimum((0, 0, 1.), axis=(1, 2))
    ion_pos[np.abs(ion_pos)<1e-6] = 0
    print('Null point found at: ', ion_pos, 'µm')
    ion_height = ion_pos[2]

    derivs = get_electrode_coeffs(s=s, ion_pos=ion_pos, filename=os.path.join(rundir, coeff_filename))
    '''
    Load saved coefficients from file and prepare to solve voltages.
    '''
    fc = load_coeffs(filename=os.path.join(rundir, coeff_filename))
    '''
    coeff_indices:
    indices of which elements from target curvature to use during fit.
    '''

    '''
    Solve for basis axial set.
    '''
    coeff_indices = np.arange(5)
    electrode_v, group_v = solve_voltages(el_names=s.names,
                    fitted_coeffs=fc, target_coeffs=target_axial_coeffs,
                    coeff_indices=coeff_indices,
                    groups=electrode_config.get('axial_groups'),
                    filename=os.path.join(rundir, Vs_axial_filename),
                    save_file=True)
    plot_length = 50.
    plot_fitted_coeffs(s=s, electrode_voltages=electrode_v,
                        target_coeffs=target_axial_coeffs,
                        ion_height=ion_height, length=plot_length,
                        shift={'z': ion_height},
                        plot_target=True)
    '''
    Solve for basis tilt set.
    '''
    r_tilt = get_axes_unitv_potentials(s=s, tilt_theta=tilt_angle_radians,
                                       shift={'z': ion_height})

    '''
        Solve for tilt curvatures due to unit voltage on each electrode.
        Use polynomial fitting to extract up to 2nd order coefficients
        from potential expansion.
    '''
    tilt_fc, residuals = get_electrode_coeffs_fit(s, *r_tilt,
                                ion_height=ion_height,
                                filename=os.path.join(rundir, tilt_coeff_filename), plot=False)
    tilt_fc = load_coeffs(filename=os.path.join(rundir, tilt_coeff_filename))
    tilt_coeff_indices = np.arange(6)
    tilt_electrode_v, tilt_group_v = solve_voltages(el_names=s.names,
                    fitted_coeffs=tilt_fc, target_coeffs=target_tilt_coeffs,
                    coeff_indices=tilt_coeff_indices,
                    groups=electrode_config.get('tilt_groups'),
                    filename=os.path.join(rundir, Vs_tilt_filename),
                    exact=False, save_file=True)
    plot_length = 20.
    plot_fitted_coeffs(s=s, electrode_voltages=tilt_electrode_v,
                        target_coeffs=target_tilt_coeffs,
                        ion_height=ion_height, length=plot_length,
                        shift={'z': ion_height},
                        plot_target=False)

    solve_freqs(s=s, f_rad=[f_rad_y, f_rad_z],
                f_axial=f_axial,
                f_traprf=f_traprf,
                m=m, q=q, l=l,
                dc_basis_axial_ref_coeffs=target_axial_coeffs,
                dc_basis_tilt_ref_coeffs=target_tilt_coeffs,
                dc_basis_axial_set_file=os.path.join(rundir, Vs_axial_filename),
                dc_basis_tilt_set_file=os.path.join(rundir, Vs_tilt_filename),

                dc_scaledfreq_axial_set_file=os.path.join(rundir, Vs_scaledfreq_axial_filename),
                dc_scaledfreq_tilt_set_file=os.path.join(rundir, Vs_scaledfreq_tilt_filename),
                dc_scaledfreq_overall_set_file=os.path.join(rundir, Vs_scaledfreq_overall_filename),

                dc_artiq_axial_set_file=os.path.join(rundir, Vs_artiq_axial_filename),
                dc_artiq_tilt_set_file=os.path.join(rundir, Vs_artiq_tilt_filename),

                axial_electrodes=axial_electrodes,
                artiq_config=artiq_config,
                translate_well=(True, translate_well_by),
                make_plot=True,
                save_result=True)

    '''
    Show data in GUI.
    '''
    ion_height_label.configure(text=f'Ion position above trap surface: {round(ion_height, 2)} µm')

row += 1
button_simulate = ttk.Button(window, 
                        text = "Simulate",
                        command = run_simulation,
                        default='active')
button_simulate.grid(column=1, row=row, **options)

ion_height_label = Label(window, 
                text = 'Ion position above trap surface: --.-- µm',
                height = 1)
ion_height_label.grid(column=2, row=row, **options)

# Let the window wait for any events
window.mainloop()