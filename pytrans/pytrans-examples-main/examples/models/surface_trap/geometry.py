#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 04/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>


# trap geometry parameters
um = 1e-6
rf_width = 120 * um
rf_sep = 60 * um

n_dc_lines = 5
dc_width = 100 * um

dc_length = 1300 * um
trap_length = 3000 * um

filename = "surface_trap_geometry.json"

corners = {}

corners['DCintop'] = (-trap_length / 2, trap_length / 2, 0, rf_sep / 2)
corners['DCinbot'] = (-trap_length / 2, trap_length / 2, -rf_sep / 2, 0)

dc_edges = [
    (dc_width * (j - n_dc_lines / 2), dc_width * (j + 1 - n_dc_lines / 2))
    for j in range(n_dc_lines)
]

corners.update({
    f"DCtop{j + 1}": dc_edges[j] + ((rf_sep / 2 + rf_width), (rf_sep / 2 + rf_width + dc_length))
    for j in range(n_dc_lines)
})

corners.update({
    f"DCbot{j + 1}": dc_edges[j] + (-(rf_sep / 2 + rf_width + dc_length), -(rf_sep / 2 + rf_width))
    for j in range(n_dc_lines)
})


def ring_from_corners(corners):
    x1, x2, y1, y2 = map(lambda x: x / um, corners)
    ring = [[x1, y1], [x2, y1], [x2, y2], [x1, y2], [x1, y1]]
    return ring


if __name__ == '__main__':
    import json
    import matplotlib.pyplot as plt
    from matplotlib.path import Path
    from matplotlib.patches import PathPatch

    json_geometry = {name: [ring_from_corners(_corners)] for name, _corners in corners.items()}

    rf_bot_corners = (-trap_length / 2, trap_length / 2, -(rf_sep / 2 + rf_width), -rf_sep / 2)
    rf_top_corners = (-trap_length / 2, trap_length / 2, rf_sep / 2, (rf_sep / 2 + rf_width))

    json_geometry["RF"] = [ring_from_corners(rf_bot_corners) + ring_from_corners(rf_top_corners)]

    print(json_geometry.keys())
    with open(filename, 'w') as f:
        json.dump(json_geometry, f, indent=2)

    fig, ax = plt.subplots()
    for name, ele in json_geometry.items():
        path = Path.make_compound_path(*[Path(ring) for ring in ele])
        patch = PathPatch(path, linewidth=1, facecolor='none')
        ax.add_patch(patch)

    ax.relim()
    ax.autoscale_view()

    plt.show()
