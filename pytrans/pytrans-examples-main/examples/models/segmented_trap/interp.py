#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 04/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>

from tqdm import tqdm
import numpy as np
from qinterp import TricubicScalarInterpolator
from .data import data_path

data = np.load(data_path / 'segmented_trap_data.npz')

X = data["X"]
Y = data["Y"]
Z = data["Z"]

electrodes = [
    "DCCa6", "DCCa7", "DCCa8",
    "DCCc6", "DCCc7", "DCCc8",
    "RF_pondpot_1V1MHz1amu"
]


class TricubicScalarInterpolator1(TricubicScalarInterpolator):

    def _ravel_coords(self, *args):
        args = np.broadcast_arrays(*args)
        shape = args[0].shape
        args = list(map(np.ravel, args))
        X = np.stack(args, axis=1).astype(float)
        return shape, X

    def __call__(self, x, y, z, d=0):
        shape, X = self._ravel_coords(x, y, z)
        return super().__call__(X, d).reshape(shape + (d * (3,)))


interpolators = {}

for name in tqdm(electrodes, "Loading interpolators"):
    field = np.stack([X, Y, Z, data[name]], axis=1)
    interp = TricubicScalarInterpolator1(field)
    interpolators[name] = interp

x_min = interp.xIntMin
x_max = interp.xIntMax
y_min = interp.yIntMin
y_max = interp.yIntMax
z_min = interp.zIntMin
z_max = interp.zIntMax
