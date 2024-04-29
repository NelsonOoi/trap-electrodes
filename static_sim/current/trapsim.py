#########################################################################
#           Functions for electrode simulation & analysis                #
##########################################################################
import matplotlib.pyplot as plt, numpy as np, scipy.constants as ct
from matplotlib.patches import Polygon
import matplotlib.cm as cm
from electrode import (System, PolygonPixelElectrode, euler_matrix,
    PointPixelElectrode, PotentialObjective,
    PatternRangeConstraint, shaped, utils)
import scipy.optimize as sciopt
from scipy.signal import argrelextrema
from scipy.optimize import minimize, linprog
from scipy.integrate import solve_ivp
from scipy.fft import fft, fftfreq
import math
import json
import pandas as pd
from datetime import datetime, date
import gdspy
import shapely.geometry as sg
import shapely.ops as so
import shapely as sh
import time
import os

def append_filepath(filename):
    filename = str(os.path.dirname(os.path.abspath(__file__))) + '/' + filename
    return filename

''' Plot trap using GDS '''
default_electrode_mapping = {
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
default_trap_center = [(4920+5040)/2, (5502.5 + 5297.5)/2]
m_to_micron = 1e6
micron_to_m = 1e-6

electrode_ordering = [str(i) for i in range(1, 21)]
list2 = ['r', 'gnd']

alternate_ordering = False

if (alternate_ordering):
    electrode_ordering = [str(i) for i in range(2, 21)]
    electrode_ordering.remove('11')
    list2 = ['1', '11', 'r', 'gnd']

electrode_ordering = electrode_ordering + list2

'''
Loads trap from GDS.

Returns:
System of electrodes - derived from Electrode package
'''
def load_trap(filename='single_chip.gds', electrode_layer=37,
              ito_layer=12, datatype=0, plot=False,
              xlim=(1000,10000), ylim=(2420,8000),
              electrode_mapping=default_electrode_mapping,
              electrode_ordering=electrode_ordering, trap_center=default_trap_center, buildup=False):
    filename = str(os.path.dirname(os.path.abspath(__file__))) + '/' + filename
    trap = gdspy.GdsLibrary(infile=filename)
    main_cell = trap.top_level()[0]
    pol_dict = main_cell.get_polygons(by_spec=True)
    # electrode patterned layer = 37 for quetzal
    electrode_polygons = pol_dict[(electrode_layer, datatype)]
    ito_polygons = []
    el_list_prev = []
    el_list_current = []
    el_dict = {}
    if(ito_layer is not None):
        ito_polygons = pol_dict[(ito_layer, datatype)]
    
    # iterative merge.
    # treat ito and electrode parts as the same.
    shapes = [*electrode_polygons, *ito_polygons]

    # cast to shapely polygons
    for shape in shapes:
        el = sg.Polygon(shape)
        el_list_current.append(el)
    # iterative merge.
    # only merge until no further changes occur.
    while (el_list_prev != el_list_current):
        n_el = len(el_list_current)
        el_list_prev = el_list_current
        el_list_current = []
        is_merged = np.zeros(n_el)

        # merge forward.
        # i.e. if other polygons later in the list are merged
        # in an earlier step, do not merge them later on
        for i in range(0, n_el):
            if(is_merged[i] == 0):
                el = el_list_prev[i]
                for j in range(i+1, n_el):
                    if (sh.overlaps(el, el_list_prev[j]) and is_merged[j] == 0):
                        el = so.unary_union([el, el_list_prev[j]])
                        is_merged[j] = 1
                el_list_current.append(el)

    # plot shapely polygons
    if (plot):
        plt.ion()
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
        ax.set_aspect('equal', adjustable='box')
    
    n_el = len(el_list_current)
    temp_names = np.arange(1, n_el+1)
    for i in range(n_el):
        # # el = Polygon(polygons[0], closed=False)
        # # ax.add_patch(el)
        # patch = ax.add_patch(el)
        # patch = plt.gca().add_patch(el)
        # print(sg.mapping(el)['coordinates'][0])
        # print(np.array(sg.mapping(el)['coordinates'][0]) - np.array([10000, 1000]))
        el = el_list_current[i]

        # if points defining the electrode are clockwise,
        # re-orient them so that they are CCW.
        # this is required as electrode assigns a negative area
        # to those polygons with CW vertex ordering.
        # if (el.exterior.is_ccw == False):
        el = sg.polygon.orient(el, sign=1.0)

        el_coordinates = np.array(sg.mapping(el).get('coordinates')[0]) - np.array(trap_center)
        el = sh.affinity.translate(el, xoff=-trap_center[0], yoff=-trap_center[1])

        name = str(temp_names[i])
        if (electrode_mapping != {}):
            name = str(electrode_mapping.get(name))
        
        if(name not in el_dict.keys()):
            el_dict[name] = [el_coordinates]
        else:
            el_dict.get(name).append(el_coordinates)

        if (plot):
            xs, ys = el.exterior.xy 
            ax.fill(xs, ys, alpha=0.5, fc='r', ec='none')
            d = fig.canvas.draw() 
            e = fig.canvas.flush_events()
            if (buildup):
                time.sleep(0.1)
            # ax.text(xs[0], ys[0], name, size=5)
    
    electrodes = list(el_dict.items())
    if(electrode_ordering != []):
        electrodes = []
        for name in electrode_ordering:
            electrodes.append((name, el_dict.get(name)))
    if(plot):
        plt.ioff()
        plt.show()

    s = System([PolygonPixelElectrode(name=n, paths=map(np.array, p))
                for n, p in electrodes])
    # plt.figure()
    # fig, ax = plt.subplots(1, 2, figsize=(15,10))
    # s.plot(ax[0])
    # plt.show()
    return s, electrodes, el_dict

spacing = 5.
n = 9
dc_width = 120.
# dc_height = 1000.
gnd_width = 300.
gnd_height = 1000.
dc_height = 1000.
dc_mid_height = 20.
# dc_mid_width = n * dc_width + (n-1) * spacing
dc_mid_width = n * dc_width + (n+1) * spacing + 2 * gnd_width
rf_height = 70.
# rf_width = n * dc_width + (n-1) * spacing
rf_width = n * dc_width + (n+1) * spacing + 2 * gnd_width

leftmost_electrode_x = -(4.5*dc_width + 4*spacing)
rightmost_electrode_x = (4.5*dc_width + 4*spacing)

gnd_start_pos = [(leftmost_electrode_x-gnd_width-spacing, (2.5*spacing + dc_mid_height + rf_height)),
                (leftmost_electrode_x-gnd_width-spacing, -(2.5*spacing + dc_mid_height + rf_height + dc_height)),
                (rightmost_electrode_x+spacing, (2.5*spacing + dc_mid_height + rf_height)),
                (rightmost_electrode_x+spacing, -(2.5*spacing + dc_mid_height + rf_height + dc_height))
                ]

dc_start_pos = [(leftmost_electrode_x, (2.5*spacing + dc_mid_height + rf_height)),
                (leftmost_electrode_x, -(2.5*spacing + dc_mid_height + rf_height + dc_height))]
rf_start_pos = [(leftmost_electrode_x-gnd_width-spacing, (1.5*spacing + dc_mid_height)),
                (leftmost_electrode_x-gnd_width-spacing, -(1.5*spacing + dc_mid_height + rf_height))]
dc_mid_start_pos = [(leftmost_electrode_x-gnd_width-spacing, (.5*spacing)),
                (leftmost_electrode_x-gnd_width-spacing, -(.5*spacing +dc_mid_height))]

def trap(spacing=5.,
        gnd_width=300., gnd_height=1000., gnd_start_pos=gnd_start_pos,
        dc_width=dc_width, dc_height=dc_height, dc_start_pos=dc_start_pos,
        rf_width=rf_width, rf_height=rf_height, rf_start_pos=rf_start_pos,
        dc_mid_width=dc_mid_width, dc_mid_height=dc_mid_height, dc_mid_start_pos=dc_mid_start_pos,
        n=9):
    '''
    Defines an approximate Quetzal trap.
    # n electrodes per row (10 default)
    # number of rows (2 default)
    # width of electrode
    # height of electrode
    # spacing between electrodes
    # tuple of coordinates to start each row
        # top row x-coordinate (bottom left point)
        # top row y-coordinate (bottom left point)
    # bottom row x-coordinate (bottom left point)
    # bottom row y-coordinate (bottom left point)
    '''
    electrodes = []

    #add dc electrodes
    #define effective size, without dielectric gaps
    w_eff = dc_width + spacing
    h_eff = dc_height + spacing
    #start with top row
    row = 0
    for sp in dc_start_pos:
        xp_start = sp[0] - spacing/2
        yp = sp[1] - spacing/2
        for electrode in range(n):
            # e.g. electrode 2 would be leftmost in top row
            # e.g. electrode 12 would be leftmost in bottom row
            electrode_id = str(row * (n+1) + electrode + 2)
            x_shift = electrode * w_eff
            xp = xp_start + x_shift
            # print(xp)
            electrodes.append(
                (electrode_id, [[
                    (xp, yp),
                    (xp + w_eff, yp),
                    (xp + w_eff, yp + h_eff),
                    (xp, yp + h_eff)
                ]])
            )
        row += 1
    dc_mid_w_eff = dc_mid_width + spacing
    dc_mid_h_eff = dc_mid_height + spacing
    for i in range (2):
        xp = dc_mid_start_pos[i][0] - spacing/2
        yp = dc_mid_start_pos[i][1] - spacing/2
        dc_mid_id = str(10*i + 1)
        electrodes.insert(10*i,
            (dc_mid_id, [[
                (xp, yp),
                (xp + dc_mid_w_eff, yp),
                (xp + dc_mid_w_eff, yp + dc_mid_h_eff),
                (xp, yp + dc_mid_h_eff)
            ]])
        )
    
    # add rf electrodes
    rf_w_eff = rf_width + spacing
    rf_h_eff = rf_height + spacing
    rf_electrodes = []
    for i in range (2):
        xp = rf_start_pos[i][0] - spacing/2
        yp = rf_start_pos[i][1] - spacing/2
        # rf_e_id = str('r'+i+1) #use this for unique rf electrode id
        rf_electrodes.append([
            (xp, yp),
            (xp + rf_w_eff, yp),
            (xp + rf_w_eff, yp + rf_h_eff),
            (xp, yp + rf_h_eff)
        ])
        rf_e_id = str('r')
    electrodes.append(
        (rf_e_id, rf_electrodes)
    )
    
    # add gnd electrodes
    #define effective size, without dielectric gaps
    gnd_w_eff = gnd_width + spacing
    gnd_h_eff = gnd_height + spacing
    gnd_electrodes = []
    for sp in gnd_start_pos:
        xp = sp[0] - spacing/2
        yp = sp[1] - spacing/2
        for gnd in range(n):
            gnd_id = 'gnd'
            gnd_electrodes.append([
                (xp, yp),
                (xp + gnd_w_eff, yp),
                (xp + gnd_w_eff, yp + gnd_h_eff),
                (xp, yp + gnd_h_eff)
            ])
    electrodes.append(
        (gnd_id, gnd_electrodes)
    )
    
    s = System([PolygonPixelElectrode(name=n, paths=map(np.array, p))
                for n, p in electrodes])
    # apply rf potential to electrode
    # s['r'].rf = 1.
    save_trap(electrodes)
    return s

def save_trap(el_map, name='quetzal'):
    '''
    Saves trap as a JSON file.
    Keys are electrode names.
    Values are electrode polygons.
    '''
    el_dict = dict(el_map)
    # Convert and write JSON object to file
    with open(str(name) + ".json", "w") as outfile: 
        json.dump(el_dict, outfile)

def get_electrode_coeffs(s, ion_pos,
                         save_file=True,
                         filename=str(datetime.now().strftime('%Y-%m-%d-%H-%M-%S')) + '_gds_quetzal.csv'):
    '''
    NOTE: This is an experimental implementation using electrode only
    to determine the curvatures, and by extension, coefficients.
    Fitting is NOT used.

    Define axes and retrieve potentials along each axis.

    Inputs:
    - System of electrodes
    - Ion position in 3D.
    '''
    filename = str(os.path.dirname(os.path.abspath(__file__))) + '/' + filename

    d0 = s.individual_potential(ion_pos, 0)
    d0 = d0[:, :, 0]

    # First derivative: E-field due to unit potential
    d1 = s.individual_potential(ion_pos, 1)
    d_x, d_y, d_z = d1[:, :, 0], d1[:, :, 1], d1[:, :, 2]

    # this yields the curvatures in an unexpanded form.
    # index 0: d_xx
    # index 3: d_yy
    # solve for d_zz using laplace.

    d2 = s.individual_potential(ion_pos, 2)
    d_xx, d_yy = d2[:, :, 0] / 2, d2[:, :, 3] / 2
    d_zz = -(d_xx + d_yy)

    # TODO: the entry with index 4 corresponds to dy dz,
    # does this also match the curvature along tilted axes?

    fit_coeffs = np.concatenate((d_x, d_xx, d_y, d_yy, d_z, d_zz), axis=1)

    full_coeffs = np.concatenate((np.array([s.names]).T,
                             d0, d_x, d_xx,
                             d0, d_y, d_yy,
                             d0, d_z, d_zz), axis=1)
    
    if (save_file):
        fit_order = 2
        axes = ['x', 'y', 'z']
        columns = ['id']
        for i in axes:
            for j in range(fit_order+1):
                columns.append('c_'+str(i)+str(j))
        df = pd.DataFrame(full_coeffs, columns=columns)
        df.to_csv(filename)
    return fit_coeffs

def get_axes_unitv_potentials(s, length=10., res=10001, shift={'z': 0}, tilt_theta=0):
    '''
    Define axes and retrieve potentials along each axis.

    Inputs:
    - System of electrodes
    - Length of axis
    - Resolution along each axis
    - Directional shift
    - Tilt of y,z-axes
    '''

    # create fit axes
    x3d, x = single_axis('x', bounds=(-length/2,length/2), res=res, shift=shift)
    y3d, y = single_axis('y', bounds=(-length/2,length/2), res=res, shift=shift)
    z3d, z = single_axis('z', bounds=(-length/2,length/2), res=res, shift=shift)
    # z = z + shift['z']
    if(tilt_theta != 0):
        y3d, y = tilt_single_axis(tilt_theta, 'y_prime', bounds=(-length/2,length/2), res=res, shift=shift)
        z3d, z = tilt_single_axis(tilt_theta, 'z_prime', bounds=(-length/2,length/2), res=res, shift=shift)
        print('y3d:', y3d)
        print('y:', y)
        print('z3d:', z3d)
        print('z:', z)
    # find individual potential contributions
    p_x = s.individual_potential(x3d, 0)
    p_y = s.individual_potential(y3d, 0)
    p_z = s.individual_potential(z3d, 0)

    return [x, y, z, p_x, p_y, p_z]

def get_axes_potentials(s, length=20., res=10001, shift={'z': 0}, tilt_theta=0):
    '''
    Define axes and retrieve potentials along each axis.

    Inputs:
    - System of electrodes
    - Length of axis
    - Resolution along each axis
    - Directional shift
    - Tilt of y,z-axes
    '''

    # create fit axes
    x3d, x = single_axis('x', bounds=(-length/2,length/2), res=res, shift=shift)
    y3d, y = single_axis('y', bounds=(-length/2,length/2), res=res, shift=shift)
    # z3d, z = single_axis('z', bounds=(np.max([0., -length/2]),length/2), res=res, shift=shift)
    z3d, z = single_axis('z', bounds=(-length/2,length/2), res=res, shift=shift)

    if(tilt_theta != 0):
        y3d, y = tilt_single_axis(tilt_theta, 'y_prime', bounds=(-length/2,length/2), res=res, shift=shift)
        z3d, z = tilt_single_axis(tilt_theta, 'z_prime', bounds=(-length/2,length/2), res=res, shift=shift)

    # find individual potential contributions
    p_x = s.electrical_potential(x3d, typ='dc', derivative=0)
    p_y = s.electrical_potential(y3d, typ='dc', derivative=0)
    p_z = s.electrical_potential(z3d, typ='dc', derivative=0)

    # z = z + shift['z']
    return [x, y, z, p_x, p_y, p_z]

'''
Helper function for solving polynomal coefficients.
'''
def solve_axes_coeffs(axis, p, order=0):
    '''
    Fit 2nd order polynomial for 1 electrode along 1 axis
    f = c + a_x1 * x + a_x2 * x^2
    Returns:
    - c_arr: Coefficients
    - A: Matrix of powers of x used in fitting
    - residuals: Error of fit.
    '''
    # create matrix of coordinate values
    c = np.ones(axis.shape[0])
    A = np.column_stack([c, axis])
    for i in range(2, order+1):
        A = np.column_stack([A, axis**i])
    c_arr = []
    c, residuals, rank, s = np.linalg.lstsq(A, p, rcond=None)
    c_arr.append(c)
    return c_arr, A, residuals[0]

def get_electrode_coeffs_fit(s, x, y, z, p_x, p_y, p_z, ion_height,
                             filename=str(datetime.now().strftime('%Y-%m-%d-%H-%M-%S')) + '_gds_quetzal.csv',
                             plot=False):
    fitted_coeffs = []
    residuals = []
    z_new = z
    axes = {'x': x, 'y': y, 'z': z}
    fit_order = 2
    for electrode_num in range(len(s.names)):
        electrode_coeffs = []
        electrode_residuals = []
        # measured = {'x': p_x[electrode_num].flatten(), 'y': p_y[electrode_num].flatten()}
        measured = {'x': p_x[electrode_num].flatten(), 'y': p_y[electrode_num].flatten(), 'z': p_z[electrode_num].flatten()}
        num_axes = len(axes)
        if (plot):
            fig_width = 23
            fig, ax = plt.subplots(1, num_axes, figsize=(fig_width, 1/num_axes*fig_width))
            ax = ax.flat
            ax[0].set_ylabel('V')
            electrode_name = s.names[electrode_num]
            fig.suptitle(f'Electrode {electrode_name}')
        for i in range(num_axes):
            ax_name = list(axes.keys())[i]
            ax_coord = list(axes.values())[i]
            # fit
            c_arr, A, residual = solve_axes_coeffs(ax_coord, measured[ax_name], order=fit_order)
            residual_str = np.format_float_scientific(residual, unique=False, precision=4)
            electrode_coeffs.append(c_arr)
            electrode_residuals.append(residual)
            # get predicted result
            predicted = np.dot(A, np.array(c_arr).T)
            if (plot):
                # plot numerical solution
                ax[i].set_xlabel(f'{ax_name} (µm)')
                if(i == 2):
                    ax[i].set_title(f'Potential along {ax_name}')
                else:
                    ax[i].set_title(f'Potential along {ax_name}, at z = {ion_height}µm')
                ax[i].plot(ax_coord, measured[ax_name], label='Measured')
                # plot predicted fit
                ax[i].plot(ax_coord, predicted, label=f'Predicted, residual: {residual_str}')
                ax[i].grid()
                ax[i].legend()
        if(plot):
            plt.show()
        fitted_coeffs.append(np.array(electrode_coeffs).flatten())
        residuals.append(np.array(electrode_residuals))

    # generate columns
    columns = ['id']
    for i in axes:
        for j in range(fit_order+1):
            columns.append('c_'+str(i)+str(j))
    for i in axes:
        columns.append('residuals_' + str(i))
    # generate pandas csv file
    # print(overall_residuals)
    savefile = np.concatenate((np.array([s.names]).T, np.array(fitted_coeffs), np.array(residuals)), axis=1)
    df = pd.DataFrame(savefile, columns=columns)
    df.to_csv(append_filepath(filename))
    # print(columns)
    return [fitted_coeffs, residuals]

def load_coeffs(filename):
    # filename = str(os.path.dirname(os.path.abspath(__file__))) + '/' + filename
    filename = append_filepath(filename)
    df = pd.read_csv(filename, usecols=['c_x1', 'c_x2', 'c_y1', 'c_y2', 'c_z1', 'c_z2'])
    return df.to_numpy()

def solve_voltages(el_names, fitted_coeffs, target_coeffs, groups, filename,
                   coeff_indices=np.arange(6), save_file=True, exact=False):
    '''
    Solves for electrode potentials.
    Inputs:
    - List of electrode names.
    - Fitted coefficients under unit voltage electrodes.
    '''
    filename = str(os.path.dirname(os.path.abspath(__file__))) + '/' + filename
    A = np.array(fitted_coeffs)
    # remove constant terms from matrix
    # A = np.delete(A, [0, 3, 6], 1)
    # remove rf and gnd electrodes
    # A = np.delete(A, [-1, -2], 0)
    # convert from 2nd order target curvature
    # to target coefficients

    fit_coeffs = []
    A_grouped = []

    for index in coeff_indices:
        fit_coeffs.append(target_coeffs[index])

    fit_coeffs = np.array(fit_coeffs)
    print('Target fit coefficients: ', fit_coeffs)
    print()

    for group in groups:
        group_coeffs = np.zeros(6)
        for electrode_name in group:
            group_coeffs += np.array(A[el_names.index(electrode_name)])

        do_rounding = False
        if (do_rounding):
            # checking if we should remove first-order
            f_order = [0, 2, 4]
            for f in f_order:
                # if (np.abs(group_coeffs[f] * m_to_micron) < threshold):
                group_coeffs[f] = np.round(group_coeffs[f] * m_to_micron, 3) / m_to_micron
                # group_coeffs[f] = np.round(group_coeffs[f] * m_to_micron, 3)
            # checking if we should remove second-order
            s_order = [1, 3, 5]
            for s in s_order:
                # if (np.abs(group_coeffs[s] * m_to_micron**2) < threshold):
                #     group_coeffs[s] = 0
                group_coeffs[f] = np.round(group_coeffs[f] * m_to_micron**2, 3) / (m_to_micron**2)
                # group_coeffs[s] = np.round(group_coeffs[s] * m_to_micron**2, 3)

        # do_clean_data = False
        # if (do_clean_data):
        #     threshold = 1e-4
        #     # checking if we should remove first-order
        #     f_order = [0, 2, 4]
        #     for f in f_order:
        #         if (np.abs(group_coeffs[f] * m_to_micron) < threshold):
        #             group_coeffs[f] = 0
        #     # checking if we should remove second-order
        #     s_order = [1, 3, 5]
        #     for s in s_order:
        #         if (np.abs(group_coeffs[s] * m_to_micron**2) < threshold):
        #             group_coeffs[s] = 0

        A_grouped.append(group_coeffs)
    
    A_grouped = np.array(A_grouped)
    '''
    A_grouped:
    An n x 6 matrix where n is the number of groups.
    '''

    ''' adjust target coefficients. '''
    # f_order = [0, 2, 4]
    # for f in f_order:
    #     target_coeffs[f] *= m_to_micron
    # s_order = [1, 3, 5]
    # for s in s_order:
    #     target_coeffs[s] *= m_to_micron ** 2
    '''
    Reduces the number of degrees of freedom.
    '''
    A_grouped_compact = np.array([A_grouped[:, coeff_indices[0]]])

    for i in range (1, len(coeff_indices)):
        A_grouped_compact = np.append(A_grouped_compact, np.array([A_grouped[:, coeff_indices[i]]]),
                                           axis=0)
    # print('A grouped compact:')
    # print(A_grouped_compact)
    # print()
    group_voltages = []
    if (exact):
        group_voltages = np.linalg.solve(A_grouped_compact, target_coeffs)
    else:
        # group_voltages, residuals, rank, s = np.linalg.lstsq(A_grouped_compact, fit_coeffs, rcond=-1)

        '''
        Test using scipy's bounded least squares,
        to keep electrode voltages within -10V < v < +10V.
        '''
        res = sciopt.lsq_linear(A=A_grouped_compact, b=fit_coeffs, bounds=(-4, 4), lsmr_tol='auto', verbose=1)
        group_voltages = res.x
    print('Matrix multiplication results: ', np.matmul(A_grouped_compact, group_voltages))

    electrode_voltages = np.zeros(len(el_names))
    for i in range(len(groups)):
        for electrode_name in groups[i]:
            electrode_voltages[el_names.index(electrode_name)] = group_voltages[i]

    # save generated voltages to csv
    if (save_file):
        el_df = pd.DataFrame(electrode_voltages.T, columns=["V"])
        el_df.to_csv(filename, index=False)
    return electrode_voltages, group_voltages

def plot_fitted_coeffs(s, electrode_voltages, target_coeffs, ion_height,
                          length=20, res=10001, shift={'z': 0}, target_coeff_indices=[1, 3, 5],
                          plot_target=True):
    x3d, x = single_axis('x', bounds=(-length/2,length/2), res=res, shift=shift)
    y3d, y = single_axis('y', bounds=(-length/2,length/2), res=res, shift=shift)
    z3d, z = single_axis('z', bounds=(-length/2,length/2), res=res, shift=shift)
    with s.with_voltages(dcs=electrode_voltages):
        vx = s.electrical_potential(x3d).T
        vy = s.electrical_potential(y3d).T
        vz = s.electrical_potential(z3d).T
        
        fig, ax = plt.subplots(1, 3, figsize=(23,7))
        ax = ax.flat
        ax[0].plot(x, vx[0], label='Electrode contribution')
        ax[0].set_title('x-axis')
        ax[1].plot(y, vy[0], label='Electrode contribution')
        ax[1].set_title('y-axis')
        ax[2].plot(z+ion_height, vz[0], label='Electrode contribution')
        ax[2].set_title('z-axis')
        ax[2].axvline(x=ion_height, linestyle='--', color='r')
        if (plot_target):
            ax[0].plot(x, target_coeffs[target_coeff_indices[0]] * x**2 +vx[0][int((vx.size+1)/2)], label='Allcock curvature')
            ax[1].plot(y, target_coeffs[target_coeff_indices[1]] * y**2 +vy[0][int((vy.size+1)/2)], label='Allcock curvature')
            ax[2].plot(z+ion_height, target_coeffs[target_coeff_indices[2]] * z**2 +vz[0][int((vz.size+1)/2)], label='Allcock curvature')
        ax[0].set_ylabel('V')
        for i in range(3):
            ax[i].set_xlabel('µm')
            ax[i].legend()
            ax[i].grid()

def solve_freqs(s, f_rad=3e6, f_split=0., f_axial=1e6, f_traprf=30e6,
                m=40*ct.atomic_mass, q=1*ct.elementary_charge,
                l=1e-6, dc_axial_ref_coeffs=[0., 2e-6, 0., -1e-6, 0., -1e-6],
                dc_tilt_ref_coeffs=[0., 0., 0., -2e-6, 0., 2e-6],
                dc_axial_set_file="Vs_2024-02-18_axial.csv",
                dc_tilt_set_file="Vs_2024-02-18_tilt.csv",
                do_plot_potential=True,
                save_result=True):

    '''
    Solves for frequencies in axial and radial directions.
    '''

    # dc_axial_set_file = str(os.path.dirname(os.path.abspath(__file__))) + '/' + dc_axial_set_file
    # dc_tilt_set_file = str(os.path.dirname(os.path.abspath(__file__))) + '/' + dc_tilt_set_file
    # overall_dc_set_file = str(os.path.dirname(os.path.abspath(__file__))) + '/'

    # 3 MHz radial frequency default target
    o_rad = 2 * np.pi * f_rad
    # 1 MHz axial frequency default target
    # 30MHz trap rf default
    o_traprf = 2 * np.pi * f_traprf
    # 40Ca+ ion mass
    # +1 ion charge
    # 1µm length scale

    ### 1. solve the DC curvature needed for 1MHz
    o_axial = 2 * np.pi * f_axial
    u_axial = o_axial**2 * m * l**2 / q
    print("Target axial frequency:", f_axial, "Hz")
    print("DC axial curvature needed", u_axial, "V/µm^2")

    # by gauss's law, laplacian of dc potential
    # in free space is 0. use (+2, -1, -1) distribution.
    # have solved for the dc voltage sets previously.
    radial_anticonfinement = u_axial / 2
    # dc curvature scales linearly with voltage.
    dc_div = u_axial/(2*np.array(dc_axial_ref_coeffs)[1])
    v_axial_scaling = np.max(dc_div[np.isfinite(dc_div)])
    print("DC axial voltage scales up by:", v_axial_scaling)
    print("\n")

    # read axial basis voltage set
    dc_axial_set, dc_tilt_set = read_electrode_voltages([dc_axial_set_file, dc_tilt_set_file])
    # scale up by required dc scaling
    print("Basis axial DC voltage set:\n", dc_axial_set.T)
    dc_axial_set *= v_axial_scaling
    print("Scaled axial DC voltage set:\n", dc_axial_set.T)

    ### 2. solve required pp curvature for 3MHz
    # with additional anticonfinement
    # and corresponding v_peak

    # assume (ux, uy, uz) is the pp curvature needed for 3MHz
    # radial frequencies. need to scale to:
    # (ux', uy+dc_anticonf, uz+dc_anticonf)

    u_rad = o_rad**2 * m * l**2 / q
    print("Target radial frequency:", f_rad, "Hz")
    print("Radial curvature needed:", u_rad, "V/µm^2")

    u_rad += radial_anticonfinement
    print("Curvature with anticonfinement compensation needed:", u_rad, "V/µm^2")
    print("\n")

    # find curvature produced at v_p = 10V

    vp = 10.
    s["r"].rf = vp*np.sqrt(q/m)/(2*l*o_traprf)
    x_guess = [0., 0., 50.]
    x0 = s.minimum(x_guess)
    print("PP saddle point found at:", x0, "µm")

    curves, modes_pp = s.modes(x0)
    curves[curves < 0] = 0
    f_rad_initial = np.sqrt(q*curves/m)/(2*np.pi*l)

    print(f"V_peak = {vp}V.")
    print("Curvature:", curves)
    print("Initial radial frequencies:", f_rad_initial)

    # since curvature scales as V_peak^2,
    # want y,z-axis curvatures to scale to target u_rad.
    # calculate using mean of y,z curvatures
    curvature_scaling = u_rad / np.mean([curves[1], curves[2]])
    print("Curvature needs scaling up by:", curvature_scaling)
    vp_scaling = np.sqrt(curvature_scaling)
    print("RF Voltage needs scaling up by:", vp_scaling)

    print("\n")
    print("Post-scaling")

    vp *= vp_scaling
    s["r"].rf = vp*np.sqrt(q/m)/(2*l*o_traprf)
    x_guess = [0., 0., 50.]
    x0 = s.minimum(x_guess)
    print("PP saddle point found at:", x0, "µm")

    curves, modes_pp = s.modes(x0)
    curves[curves < 0] = 0
    achievable_f_rad = np.sqrt(q*curves/m)/(2*np.pi*l)

    print(f"Scaled V_peak = {vp}V.")
    print("Curvature:", curves)
    print("Scaled radial frequencies:", achievable_f_rad)
    print("\n")

    # scale up by required dc scaling
    v_tilt_scaling, tilt_coeff_needed = solve_tilt_scaling(f_rad=f_rad,
                        f_split=f_split, m=m, q=q, l=l,
                        dc_tilt_ref_coeffs=dc_tilt_ref_coeffs)
    print('Tilt coeff needed:', tilt_coeff_needed)
    ### read tilt basis voltage set
    # dc_tilt_df = pd.read_csv(dc_tilt_set_file)
    # dc_tilt_set = dc_tilt_df["V"]
    # dc_tilt_set = read_electrode_voltages([dc_tilt_set_file])
    print("DC tilt voltage scales up by:", v_tilt_scaling)
    print("Basis tilt DC voltage set:\n", dc_tilt_set.T)
    dc_tilt_set *= v_tilt_scaling
    print("Scaled tilt DC voltage set:\n", dc_tilt_set.T)

    dc_electrode_set = dc_axial_set + dc_tilt_set
    print('Overall electrode set:', dc_electrode_set)
    # dc_electrode_set = dc_axial_set
    print(dc_electrode_set)
    with s.with_voltages(dcs=dc_electrode_set):
        x_guess = [0., 0., 50.]
        x0 = s.minimum(x_guess)
        ion_height = x0[2]
        tmp_len = 20.
        res = 10001
        shift = {'z': x0[2]}
        x3d, x = single_axis('x', (-tmp_len/2,tmp_len/2), res, shift=shift)
        y3d, y = single_axis('y', (-tmp_len/2,tmp_len/2), res, shift=shift)
        z3d, z = single_axis('z', (-tmp_len/2,tmp_len/2), res, shift=shift)
        vx = s.potential(x3d, 0).T
        vy = s.potential(y3d, 0).T
        vz = s.potential(z3d, 0).T
        fig, ax = plt.subplots(1, 3, figsize=(15,5))
        ax = ax.flat
        curves, modes_pp = s.modes(x0)
        print("Curvature:", curves)
        # make plots of overall rf + dc potential
        ax[0].plot(x, vx, label='Electrode contribution')
        ax[1].plot(y, vy, label='Electrode contribution')
        ax[2].plot(z + ion_height, vz, label='Electrode contribution')
        ax[2].axvline(x=ion_height, linestyle='--', color='r')
        axes_names = ['x', 'y', 'z']
        ax[0].set_ylabel('V')
        for i in range(len(ax)):
            ax[i].set_title(f'{axes_names[i]}-axis RF + DC potential')
            ax[i].set_xlabel('µm')
            ax[i].grid()

        # uncomment for static & mathieu analysis
        rf_scale = s.rf_scale(m, q, l, o_traprf)
        mu, b = s.mathieu(x0, scale=rf_scale, r=4, sorted=True)
        freqs = mu[:3].imag*o_traprf/(2*np.pi)
        modes = b[len(b)//2 - 3:len(b)//2, :3].real
        for i in range(len(freqs)):
            print(f"Mathieu frequency {i+1}:", freqs[i]/1e6, "MHz", modes[i])
        '''
        for line in s.analyze_static(x0, axis=(1, 2), m=m, q=q, l=l, o=o_traprf):
            print(line)
        '''
        if(do_plot_potential):
            plot_potential(s=s, freqs=freqs, modes=modes)
            plot_field(s=s)
            plt.show()
    
    if (save_result):
        overall_dc_df = pd.DataFrame(dc_electrode_set)
        overall_dc_set_file = f'RF{round(f_traprf/1e6, 0)}MHz-trapx{round(freqs[0]/1e6, 0)}MHz-y_prime{round(freqs[1]/1e6, 0)}MHz-z_prime{round(freqs[2]/1e6, 0)}MHz-{datetime.now().strftime("%Y-%m-%d-%H-%M-%S")}.csv'
        overall_dc_df.to_csv(append_filepath(overall_dc_set_file), index=False)

def solve_tilt_scaling(f_rad=3e6, f_split=0.2e6,
                m=40*ct.atomic_mass, q=1*ct.elementary_charge,
                l=1e-6, dc_tilt_ref_coeffs=[0., 0., 0., -1e-6, 0., +1e-6]):
    o_rad_minor = 2*np.pi*(f_rad + f_split/2) # minor axis of elliptic potential -> higher curvature
    o_rad_major = 2*np.pi*(f_rad - f_split/2) # major axis of elliptic potential
    u_rad_minor = o_rad_minor**2 * m * l**2 / q
    u_rad_major = o_rad_major**2 * m * l**2 / q
    u_diff = u_rad_minor - u_rad_major
    v_tilt_scaling = (u_diff / 2) / (2 * dc_tilt_ref_coeffs[-1])
    return v_tilt_scaling, u_diff/2/2


def plot_potential(s, d_r=10, contour_res=1000, freqs=[], modes=[]):
    '''
    Plots the potential distribution in the xy and yx planes.
    Produces contour plots from which the curvature can be intuited.
    '''
    ion_height = s.minimum([0., 0., 50.])[2]
    # ion_height = 51.7µm for GDS
    # make grid for potential view in z = z0 plane - top view looking down
    grid_xy, tsx0, tsy0 = make_xy_grid_flat([-40.,40.], [-40., 40.], ion_height, [101, 101])
    # make grid for potential view in x = 0 plane - side view
    grid_yz, tsy1, tsz1 = make_yz_grid_flat(0., [-d_r,d_r], [np.round(ion_height)-d_r, np.round(ion_height)+d_r], [101, 101])
    grids = [[grid_xy, tsx0, tsy0], [grid_yz, tsy1, tsz1]]
    
    fig, ax = plt.subplots(1, 2, figsize=(23,7))

    titles = [f'Potential in xy-plane at z={ion_height}', f'Potential in yz-plane at x={0}']
    labels = [['x-axis(µm)', 'y-axis(µm)'], ['y-axis(µm)', 'z-axis(µm)']]
    for i in range(2):
        gridpot = s.potential(grids[i][0], derivative=0)
        gridpot = np.reshape(gridpot, (len(grids[i][2]), len(grids[i][1])))
        bgc = ax[i].contourf(grids[i][1], grids[i][2], gridpot, levels=contour_res, cmap=cm.plasma)
        cp = ax[i].contour(grids[i][1], grids[i][2], gridpot, levels=10, colors='white')
        ax[i].clabel(cp, fontsize=10, colors='white')
        ax[i].set_title(titles[i])
        ax[i].set_xlabel(labels[i][0])
        ax[i].set_ylabel(labels[i][1])
        ax[i].set_aspect('equal')
        fig.colorbar(bgc)
    scaling = 10
    axial = modes[0] * scaling
    radial_1 = modes[1] * scaling
    radial_2 = modes[2] * scaling
    
    arrow_style = {
        "head_width": 1,
        "head_length": 2,
    }
    color = 'orange'
    ax[0].arrow(0., 0., dx=axial[0], dy=axial[1], color=color, **arrow_style)
    # ax[0].arrow(0., 0., dx=radial_1[0], dy=radial_1[1], color=color, **arrow_style)
    ax[0].annotate(f'{round(freqs[0]/1e6, 3)} MHz', xy=[axial[0]/2, axial[1]/2], xytext=(10, -10), color=color, textcoords='offset points')
    # ax[0].annotate(f"tilted-z': {round(freqs[1]/1e6, 3)} MHz", xy=[radial_1[0]/2, radial_1[1]/2], xytext=(10, 5), color=color, textcoords='offset points')
    
    yz_scaling = freqs[1] / (freqs[1] + freqs[2])
    yz_arrow_style = {
        "head_width": 0.5,
        "head_length": 1,
    }
    ax[1].arrow(0., ion_height, dx=radial_1[1]*yz_scaling, dy=radial_1[2]*yz_scaling, color=color, **yz_arrow_style)
    ax[1].arrow(0., ion_height, dx=radial_2[1]*(1-yz_scaling), dy=radial_2[2]*(1-yz_scaling), color=color, **yz_arrow_style)
    ax[1].annotate(f'{round(freqs[1]/1e6, 3)} MHz', xy=[radial_1[1]/2, radial_1[2]/2 + ion_height], xytext=(4, 0), color=color, textcoords='offset points')
    ax[1].annotate(f'{round(freqs[2]/1e6, 3)} MHz', xy=[radial_2[1]/2, radial_2[2]/2 + ion_height], xytext=(4, 0), color=color, textcoords='offset points')
    # plt.show()

def plot_field(s, grid=0, x_grid_bounds=(-100., 100.),
                        y_grid_bounds=(-100., 100.),
                        z_grid_bounds=(10., 100.),
                        grid_res=(101, 101, 101),
                        is_dc_potential=False,
                        do_potential_plot=False):

    grid_xz, X, Z = make_xz_grid_flat(x_grid_bounds=x_grid_bounds,
                                    # y_grid_bounds=(0., 0.),
                                    y = 0,
                                    z_grid_bounds=z_grid_bounds, grid_res=grid_res)

    grid_xz = np.array(grid_xz)

    grid_yz, Y, Z = make_yz_grid_flat(x = 0,
                                    y_grid_bounds=y_grid_bounds,
                                    z_grid_bounds=z_grid_bounds, grid_res=grid_res)

    grid_yz = np.array(grid_yz)

    grids = [grid_xz, grid_yz]
    '''
    Give labels as [x-axis, y-axis, title].
    '''
    labels = [['x (µm)', 'z (µm)', 'Electric field in xz-plane. y=0.'], ['y (µm)', 'z (µm)', 'Electric field in yz-plane. x=0.']]
    
    # for i in grid[:, 0]:
    #     print (i)

    # def make_xz_grid(x_grid_bounds, z_grid_bounds, grid_res):
    # ts_x = np.linspace(x_grid_bounds[0], x_grid_bounds[1], grid_res[0])
    # ts_z = np.linspace(z_grid_bounds[0], z_grid_bounds[1], grid_res[1])
    # grid = np.array([[(x,0,z) for x in ts_x] for z in ts_z])
    # print(grid[:, :, 2])

    # return grid, ts_x, ts_z
    # E = -dV/dr
    field = 0
    # if (is_dc_potential):
    #     field = - np.array(s.electrical_potential(grid, typ='dc', derivative=1))
    # else:
    #     field = - np.array(s.potential(grid, derivative=1))
    # grid = np.reshape(grid, (grid_res[0], grid_res[1], 3))
    # field = np.reshape(field, (grid_res[0], grid_res[1], 3))

    # fig, ax1 = plt.subplots(1, 1, figsize=(20, 20*np.sum(np.abs(z_grid_bounds)) / np.sum(np.abs(x_grid_bounds))))
    fig, ax = [0, 0]
    if (do_potential_plot):
        fig, ax = plt.subplots(1, 3, figsize=(13, 5))
    else:
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax = ax.flatten()
    for i in range(2):
        grid = grids[i]
        label = labels[i]
        if (is_dc_potential):
            field = - np.array(s.electrical_potential(grid, typ='dc', derivative=1))
            potential = np.array(s.electrical_potential(grid, typ='dc', derivative=0))
        else:
            field = - np.array(s.potential(grid, derivative=1))
            potential = np.array(s.potential(grid, derivative=0))

        # print(grid)
        grid = np.reshape(grid, (grid_res[1], grid_res[0], 3))
        # print(field, field.shape)
        # print(np.linalg.norm(field, axis=1))
        field = np.reshape(field, (grid_res[1], grid_res[0], 3))
        potential = np.reshape(potential, (grid_res[1], grid_res[0]))
        # print(grid, grid.shape)
        # print(field, field.shape)
        # print(grid[:, :, 0])
        # print(field[:, :, 0])
        # print(np.linalg.norm(field, axis=2))
        # print(field[:,:,0][0],field[:,:,1][0],field[:,:,2][0], field.shape)
        strm = ax[i].streamplot(grid[:, :, i], grid[:, :, 2], field[:, :, i], field[:, :, 2], density=0.5,
        # ax1.streamplot(grid[:, 0], grid[:,2], field[:,0], field[:,2], density=2,
        # ax1.streamplot(X, Z, field[:,0], field[:,2], density=2)
            color=np.linalg.norm(field, axis=2)*1e6, linewidth=1, cmap=cm.plasma, broken_streamlines=False)
        ax[i].set_xlabel(label[0])
        ax[i].set_ylabel(label[1])
        ax[i].set_title(label[2])
        cbar = fig.colorbar(strm.lines)
        cbar.set_label('V / m')
        if(i == 0 and do_potential_plot):
            contr_plot = ax[2].contourf(potential)
            contr_cbar = fig.colorbar(contr_plot)

    # plt.plot()

def read_electrode_voltages(files=[]):
    overall = []
    for filename in files:
        filename = append_filepath(filename=filename)
        el_v = pd.read_csv(filename, header=None)
        if ('V' in el_v.columns or 'V' in el_v.values):
            el_v = el_v.to_numpy().flatten()
            el_v = np.delete(el_v, np.where(el_v == 'V'))
            el_v = np.float64(el_v)
        else:
            el_v = el_v.to_numpy().flatten()
        if(len(el_v) == 20):
            el_v = np.append(el_v, np.array([0, 0]))
        overall.append(np.array(el_v))
    return np.array(overall)

# axis producers.

def single_axis(axis, bounds, res, shift):
    ts = np.linspace(bounds[0], bounds[1], res)
    zeros = np.zeros(res)
    line = ts
    line3d = []
    if (axis == 'x'):
        line3d = np.column_stack([ts.T, zeros, zeros])
    if (axis == 'y'):
        line3d = np.column_stack([zeros, ts.T, zeros])
    if (axis == 'z'):
        line3d = np.column_stack([zeros, zeros, ts.T])

    # to shift line height
    for key, val in shift.items():
        if (key == 'x'):
            line3d += np.column_stack([np.ones(res) * val, zeros, zeros])
        if (key == 'y'):
            line3d += np.column_stack([zeros, np.ones(res) * val, zeros])
        if (key == 'z'):
            line3d += np.column_stack([zeros, zeros, np.ones(res) * val])
    return line3d, line

def tilt_single_axis(theta, axis, bounds, res, shift={'z': 0}):
    '''
    theta in radians.
    '''
    ts = np.linspace(bounds[0], bounds[1], res)
    ts_cos = ts * np.cos(theta)
    ts_sin = ts * np.sin(theta)
    zeros = np.zeros(res)
    line = ts
    line3d = []
    if (axis == 'y_prime'):
        line3d = np.column_stack([zeros, ts_cos.T, ts_sin.T])
    if (axis == 'z_prime'):
        line3d = np.column_stack([zeros, -ts_sin.T, ts_cos.T])
        # line3d = line3d[::-1]
    # to shift line height
    for key, val in shift.items():
        if (key == 'x'):
            line3d += np.column_stack([np.ones(res) * val, zeros, zeros])
        if (key == 'y'):
            line3d += np.column_stack([zeros, np.ones(res) * val, zeros])
        if (key == 'z'):
            line3d += np.column_stack([zeros, zeros, np.ones(res) * val])
    return line3d, line

def make_xy_grid_flat(x_grid_bounds, y_grid_bounds, z, grid_res):
    ts_x = np.linspace(x_grid_bounds[0], x_grid_bounds[1], grid_res[0])
    ts_y = np.linspace(y_grid_bounds[0], y_grid_bounds[1], grid_res[1])
    grid = []
    for y in ts_y:
        for x in ts_x:
            grid.append([x, y, z])
    return grid, ts_x, ts_y

def make_yz_grid_flat(x, y_grid_bounds, z_grid_bounds, grid_res):
    ts_y = np.linspace(y_grid_bounds[0], y_grid_bounds[1], grid_res[0])
    ts_z = np.linspace(z_grid_bounds[0], z_grid_bounds[1], grid_res[1])
    grid = []
    for z in ts_z:
        for y in ts_y:
            grid.append([x, y, z])
    return grid, ts_y, ts_z

def make_xz_grid_flat(x_grid_bounds, y, z_grid_bounds, grid_res):
    ts_x = np.linspace(x_grid_bounds[0], x_grid_bounds[1], grid_res[0])
    ts_z = np.linspace(z_grid_bounds[0], z_grid_bounds[1], grid_res[1])
    grid = []
    for z in ts_z:
        for x in ts_x:
            grid.append([x, y, z])
    return grid, ts_x, ts_z

def make_xyz_grid_flat(x_grid_bounds, y_grid_bounds, z_grid_bounds, grid_res):
    ts_x = np.linspace(x_grid_bounds[0], x_grid_bounds[1], grid_res[0])
    ts_y = np.linspace(y_grid_bounds[0], y_grid_bounds[1], grid_res[1])
    ts_z = np.linspace(z_grid_bounds[0], z_grid_bounds[1], grid_res[2])
    grid = []
    for z in ts_z:
        for y in ts_y:
            for x in ts_x:
                grid.append([x, y, z])
    return grid, ts_x, ts_y, ts_z


def trapf_minima():
    # Plot - pp-frequency against rf frequency, at different u
    u_range = np.arange(2, 10+1)
    f_range = np.linspace(25e6, 50e6, 1000)
    fig, ax = plt.subplots(1, 5, figsize=(25,5))
    for u in u_range:
        l = 1e-6 # µm length scale
        m = 40*ct.atomic_mass # 40Ca+ ion mass
        q = 1*ct.elementary_charge # ion charge
        o_range = 2*np.pi*f_range
        x = [0., 0., 50.]
        pf_y = []
        pf_z = []
        x0s = []
        for o in o_range:
            s["r"].rf = u*np.sqrt(q/m)/(2*l*o)
            x0 = s.minimum(x)
            curves, modes_pp = s.modes(x0)
            freqs_pp = np.sqrt(q*curves/m)/(2*np.pi*l)
            pf_y.append(freqs_pp[1])
            pf_z.append(freqs_pp[2])
            x0s.append(x0)
            # print(x0)
            # print(curves, 'modes: ', modes_pp)
        
        ax[0].plot(f_range/1e6, pf_y, label=f'Vpeak = {u}V')
        ax[0].set_title('y-axis trap frequency (pp only)')
        # ax[0].set_xlabel('f (Hz)')
        # ax[0].set_ylabel('Curvature')
        ax[1].plot(f_range/1e6, pf_z, label=f'Vpeak = {u}V')
        ax[1].set_title('z-axis trap frequency (pp only)')

        coord = ['x', 'y', 'z']
        x0s = np.array(x0s)
        # print(x0s)
        for i in range(2, 5):
            ax[i].plot(f_range/1e6, x0s[:, i-2], label=f'Vpeak = {u}V')
            ax[i].set_title(f'Position of minima on {coord[i-2]} axis')
    
    for i in range(2):
        ax[i].set_xlabel('RF (MHz)')
        ax[i].set_ylabel('Frequency')
        ax[i].legend()
        # ax[i].set_yscale('log')
        ax[i].grid()
    for i in range(2, 5):
        ax[i].legend()
        ax[i].grid()
    for i in range(2, 4):
        ax[i].set_ylim(-0.05, 0.05)
    # plt.suptitle('Curvature against frequency due to RF only.')