from trapsim import *

s, electrodes, electrodes_dict = load_trap(filename='single_chip.gds')

z_height = 0.001
check_coords=[
    [0., 12.5, z_height],
]

main_electrode_y_coord = 2.5 + 20. + 5. + 70. + 5. + 5.
main_electrode_x_coord_start = -4. * (120. + 5.)

for i in range(9):
    check_coords.append([main_electrode_x_coord_start + i * (120. + 5.), main_electrode_y_coord, z_height])

check_coords.append([0., -12.5, z_height])

for i in range(9):
    check_coords.append([main_electrode_x_coord_start + i * (120. + 5.), -main_electrode_y_coord, z_height])

voltages = np.zeros(22)
for i in range (20):
    voltages[i-1] = 0
    voltages[i] = 5.123
    with s.with_voltages(dcs=voltages):
        pots = s.electrical_potential(check_coords, typ='dc', derivative=0)
        pots = pots.T.flatten()
        print(voltages, pots)
        print(voltages[i], pots[i])

# grid, ts_x, ts_y = make_xy_grid_flat()


