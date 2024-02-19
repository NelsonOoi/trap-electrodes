#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 04/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>

import numpy as np
from pytrans.electrode import DCElectrode, RFElectrode
from pytrans.abstract_model import AbstractTrapModel
from .interp import TricubicScalarInterpolator1, interpolators
from . import interp

import logging

logger = logging.getLogger(__name__)

_DC_ELECTRODES = [
    "DCCa6", "DCCa7", "DCCa8",
    "DCCc6", "DCCc7", "DCCc8"
]


class SegmentedTrapDCElectrode(DCElectrode):

    def __init__(self, interp_name):
        self._interp: TricubicScalarInterpolator1 = interpolators[interp_name]
        super().__init__()

    def _unit_potential(self, x, y, z):
        return self._interp(x, y, z, d=0)

    def _unit_gradient(self, x, y, z):
        return self._interp(x, y, z, d=1)

    def _unit_hessian(self, x, y, z):
        return self._interp(x, y, z, d=2)


class SegmentedTrapRFElectrode(RFElectrode):

    def __init__(self):
        self._interp: TricubicScalarInterpolator1 = interpolators['RF_pondpot_1V1MHz1amu']
        super().__init__()

    def _unit_potential(self, x, y, z):
        return self._interp(x, y, z, d=0)

    def _unit_gradient(self, x, y, z):
        return self._interp(x, y, z, d=1)

    def _unit_hessian(self, x, y, z):
        return self._interp(x, y, z, d=2)


class SegmentedTrap(AbstractTrapModel):

    _dc_electrodes = {name: SegmentedTrapDCElectrode(name) for name in _DC_ELECTRODES}
    _rf_electrodes = {"RF": SegmentedTrapRFElectrode()}
    _rf_voltage = 360.187
    _rf_freq_mhz = 113.733

    # Extra attributes and methodes to enrich the model
    w_ele = 320e-6
    y0 = 0.0
    z0 = 0.0
    dt = 10e-9
    rf_null_coords = (None, y0, z0)

    x_min = interp.x_min
    x_max = interp.x_max
    y_min = interp.y_min
    y_max = interp.y_max
    z_min = interp.z_min
    z_max = interp.z_max

    def __init__(self, use_electrodes=None):
        super().__init__(use_electrodes)
        self.x = np.linspace(self.x_min, self.x_max, 100)

    @classmethod
    def x_ele(cls, j):
        # center position of electrode j
        # j = 1 .. 5
        return cls.w_ele * (j - 7)
