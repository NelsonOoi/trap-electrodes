##########################################################################
#           Functions fpr electrode simulation & analysis                #
##########################################################################
import matplotlib.pyplot as plt, numpy as np, scipy.constants as ct
from matplotlib.patches import Polygon
import matplotlib.cm as cm
from electrode import (System, PolygonPixelElectrode, euler_matrix,
    PointPixelElectrode, PotentialObjective,
    PatternRangeConstraint, shaped, utils)
from scipy.signal import argrelextrema
from scipy.optimize import minimize
import math
import json
import pandas as pd
from datetime import datetime
from scipy.optimize import linprog
import gdspy
import shapely.geometry as sg
import shapely.ops as so
import shapely as sh
import time

def load_trap(filename='./single_chip.gds', electrode_layer=37, ito_layer=None, datatype=0, plot=True, xlim=(1000,10000), ylim=(2420,8000), show=True):
    trap = gdspy.GdsLibrary(infile=filename)
    main_cell = trap.top_level()[0]
    pol_dict = main_cell.get_polygons(by_spec=True)
    # electrode patterned layer = 37 for quetzal
    electrode_polygons = pol_dict[(electrode_layer, datatype)]

    ito_polygons = []
    if(ito_layer is not None):
        ito_polygons = pol_dict[(ito_layer, datatype)]
    if (plot):
        plt.ion()
        fig, ax = plt.subplots(figsize=(8, 8))
        # # el = Polygon(polygons[0], closed=False)
        # # ax.add_patch(el)
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
        ax.set_aspect('equal', adjustable='box')

        el_list = []
        # go through all electrodes and add any overlapping ITO 
        '''
        for electrode_part in electrode_polygons:
            # el = Polygon(shape, closed=False)
            print(electrode_part)
            el = sg.Polygon(electrode_part)

            # merge with ITO layer to form complete electrode
            for ito_part in ito_polygons:
                ito = sg.Polygon(ito_part)
                if (sh.overlaps(el, ito)):
                    el = so.unary_union([el, ito])
            print(el)
            el_list.append(el)
        '''
        
        # Merge all electrodes which are now connected after previous merge stage
        # this should be an iterative/recursive flood fill,
        # it should terminate only if no other changes occur
        # but we keep it simple for now and merge twice
        # assume that ito can at most merge two electrode parts
        # this is a reasonable assumption for quetzal
        '''
        el_list_final = []
        n_el = len(el_list)
        for i in range(0, n_el):
            el = el_list[i]
            for j in range(i+1, n_el):
                if (sh.overlaps(el, el_list[j])):
                    el = so.unary_union([el, el_list[j]])
            el_list_final.append(el)
        '''
        shapes = [*electrode_polygons, *ito_polygons]
        el_list_prev = []
        el_list_current = []
        # treat ito and electrode parts as the same. (test)
        for shape in shapes:
            el = sg.Polygon(shape)
            el_list_current.append(el)
        while (el_list_prev != el_list_current):
            print('unmatched')
            n_el = len(el_list_current)
            el_list_prev = el_list_current
            el_list_current = []
            is_merged = np.zeros(n_el)
            for i in range(0, n_el):
                if(is_merged[i] == 0):
                    el = el_list_prev[i]
                    for j in range(i+1, n_el):
                        # if (sh.overlaps(el, el_list_prev[j])):
                        if (sh.overlaps(el, el_list_prev[j]) and is_merged[j] == 0):
                            el = so.unary_union([el, el_list_prev[j]])
                            is_merged[j] = 1
                    el_list_current.append(el)
            
        
        # print(el_list_final)
        # for el in el_list_final:
        for el in el_list_current:
            # print(el)
            # print(sh.overlaps(el, el))
            xs, ys = el.exterior.xy 
            ax.fill(xs, ys, alpha=0.5, fc='r', ec='none')
            # patch = ax.add_patch(el)
            # patch = plt.gca().add_patch(el)
            d = fig.canvas.draw() 
            e = fig.canvas.flush_events()
            # name = input('electrode name:')
            if show:
                time.sleep(0.5)
        plt.ioff()
        plt.show()

def trap(spacing,
        gnd_width, gnd_height, gnd_start_pos,
        dc_width, dc_height, dc_start_pos,
        rf_width, rf_height, rf_start_pos,
        dc_mid_width, dc_mid_height, dc_mid_start_pos,
        n=9, r=2, real=False):
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
        electrodes.append(
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
    Saves tarp as a JSON file.
    Keys are electrode names.
    Values are electrode polygons.
    '''
    el_dict = dict(el_map)
    # Convert and write JSON object to file
    with open(str(name) + ".json", "w") as outfile: 
        json.dump(el_dict, outfile)

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

def solve_electrode_contributions(s, x, y, z, p_x, p_y, p_z, do_plot=True):
    overall_fitted_coeffs = []
    overall_residuals = []
    z_new = z
    axes = {'x': x, 'y': y, 'z': z}
    fit_order = 2
    for electrode_num in range(len(s.names)):
        electrode_coeffs = []
        electrode_residuals = []
        # measured = {'x': p_x[electrode_num].flatten(), 'y': p_y[electrode_num].flatten()}
        measured = {'x': p_x[electrode_num].flatten(), 'y': p_y[electrode_num].flatten(), 'z': p_z[electrode_num].flatten()}
        num_axes = len(axes)
        if (do_plot):
            fig_width = 15
            fig, ax = plt.subplots(1, num_axes, figsize=(fig_width, 1/2*fig_width))
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
            if (do_plot):
                # plot numerical solution
                ax[i].set_xlabel('µm')
                if(i == 2):
                    ax[i].set_title(f'Voltage along {ax_name}')
                else:
                    ax[i].set_title(f'Voltage along {ax_name}, at z = {z_h}µm')
                ax[i].plot(ax_coord, measured[ax_name], label='Measured')
                # plot predicted fit
                ax[i].plot(ax_coord, predicted, label=f'Predicted, residual: {residual_str}')
                ax[i].grid()
                ax[i].legend()
        overall_fitted_coeffs.append(np.array(electrode_coeffs).flatten())
        overall_residuals.append(np.array(electrode_residuals))

    # generate columns
    columns = ['id']
    for i in axes:
        for j in range(fit_order+1):
            columns.append('c_'+str(i)+str(j))
    for i in axes:
        columns.append('residuals_' + str(i))
    # generate pandas csv file
    # print(overall_residuals)
    savefile = np.concatenate((np.array([s.names]).T, np.array(overall_fitted_coeffs), np.array(overall_residuals)), axis=1)
    df = pd.DataFrame(savefile, columns=columns)
    df.to_csv('test_'+ str(datetime.now()) + '_coefficients_quetzal_test_small.csv')
    print(columns)

def solve_voltages(el_names, fitted_coeffs, b, groups):
    '''
    Solves for electrode potentials.
    Inputs:
    - List of electrode names.
    - Fitted coefficients under unit voltage electrodes.
    '''
    A = np.array(fitted_coeffs)
    # remove constant terms from matrix
    A = np.delete(A, [0, 3, 6], 1)
    # remove rf and gnd electrodes
    A = np.delete(A, [-1, -2], 0)
    b = np.array(b)

    A_grouped = []
    for group in groups:
        group_coeffs = np.zeros(6)
        for electrode_name in group:
            group_coeffs += np.array(A[el_names.index(electrode_name)])
        A_grouped.append(group_coeffs)
    A_grouped = np.array(A_grouped)
    A_grouped[np.abs(A_grouped)<1e-17] = 0

    group_voltages = np.linalg.solve(A_grouped.T, b)
    return group_voltages

def solve_freqs(s, f_rad = 3e6, f_axial = 1e6, f_traprf = 30e6,
                m = 40*ct.atomic_mass, q = 1*ct.elementary_charge,
                l = 1e-6, u_dc_ref = [0, 4e-6, 0, -2e-6, 0, -2e-6],
                dc_axial_set_file="Vs_2024-02-18_axial.csv",
                do_plot_potential=True):
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
    print("DC curvature needed", u_axial, "V/µm^2")

    # by laplace's equation, laplacian of dc potential
    # in free space is 0. use (+2, -1, -1) distribution.
    # have solved for the dc voltage sets previously.

    radial_anticonfinement = u_axial / 2

    # dc curvature scales linearly with voltage.
    u_dc = [0, u_axial, 0, -radial_anticonfinement, 0, -radial_anticonfinement]
    dc_div = np.array(u_dc)/np.array(u_dc_ref)
    vdc_scaling = np.max(dc_div[np.isfinite(dc_div)])
    print("DC Voltage needs scaling up by:", vdc_scaling)


    print("\n")

    ### read axial basis voltage set
    # scale up by required dc sclaing
    dc_df = pd.read_csv(dc_axial_set_file)
    dc_axial_set = dc_df["V"]
    print("Basis axial DC voltage set:\n", dc_axial_set.T)
    dc_axial_set *= vdc_scaling
    print("Scaled axial DC voltage set:\n", dc_axial_set.T)
    dc_axial_set = np.array(dc_axial_set)

    ### 2. solve required pp curvature for 3MHz
    # with additional anticonfinement
    # and corresponding v_peak

    # assume (ux, uy, uz) is the pp curvature needed for 3MHz
    # radial frequencies. need to scale to:
    # (ux', uy+dc_anticonf, uz+dc_anticonf)

    u_rad = o_rad**2 * m * l**2 / q
    print("Target radial frequency:", f_rad, "Hz")
    print("Curvature needed:", u_rad, "V/µm^2")

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
    f_rad = np.sqrt(q*curves/m)/(2*np.pi*l)

    print(f"Scaled V_peak = {vp}V.")
    print("Curvature:", curves)
    print("Scaled radial frequencies:", f_rad)

    print("\n")

    with s.with_voltages(dc_axial_set):
        tmp_len = 20
        res = 10000
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
        ax[2].plot(z, vz, label='Electrode contribution')

        # uncomment for static & mathieu analysis
        rf_scale = s.rf_scale(m, q, l, o_traprf)
        mu, b = s.mathieu(x0, scale=rf_scale, r=4, sorted=True)
        freqs = mu[:3].imag*o/(2*np.pi)
        modes = b[len(b)//2 - 3:len(b)//2, :3].real
        for i in range(len(freqs)):
            print(f"Mathieu frequency {i+1}:", freqs[i]/1e6, "MHz", modes[i])
        '''
        for line in s.analyze_static(x0, axis=(1, 2), m=m, q=q, l=l, o=o_traprf):
            print(line)
        '''
        if(do_plot_potential):
            plot_potential(s=s, z0=x0[2])
        print(vx)
        print(x0[2])
        print(z3d)
        
def plot_potential(s):
    '''
    Plots the potential distribution in the xy and yx planes.
    Produces contour plots from which the curvature can be intuited.
    '''
    z_h = s.minimum([0, 0, 50.])[2]
    # make grid for potential view in z = z0 plane - top view looking down
    grid_xy, tsx0, tsy0 = make_xy_grid_flat([-40.,40.], [-40., 40.], z_h, [101, 101])
    # make grid for potential view in x = 0 plane - side view
    d_r = 5.
    grid_yz, tsy1, tsz1 = make_yz_grid_flat(0., [-d_r,d_r], [np.round(z_h)-d_r, np.round(z_h)+d_r], [101, 101])
    grids = [[grid_xy, tsx0, tsy0], [grid_yz, tsy1, tsz1]]
    
    fig, ax = plt.subplots(1, 2, figsize=(20,7))
    contour_res = 1000

    titles = [f'Potential in xy-plane at z={z_h}', f'Potential in yz-plane at x={0}']
    labels = [['x-axis(µm)', 'y-axis(µm)'], ['y-axis(µm)', 'z-axis(µm)']]
    for i in range(2):
        gridpot = s.potential(grids[i][0])
        gridpot = np.reshape(gridpot, (len(grids[i][2]), len(grids[i][1])))
        bgc = ax[i].contourf(grids[i][1], grids[i][2], gridpot, levels=contour_res, cmap=cm.plasma)
        cp = ax[i].contour(grids[i][1], grids[i][2], gridpot, colors='white')
        ax[i].clabel(cp, fontsize=10, colors='white')
        ax[i].set_title(titles[i])
        ax[i].set_xlabel(labels[i][0])
        ax[i].set_ylabel(labels[i][1])
        fig.colorbar(bgc)
    plt.gca().set_aspect('equal')

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

def tilt_single_axis(theta, axis, bounds, res, shift):
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