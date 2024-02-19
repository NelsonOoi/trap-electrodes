#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Created: 01-2021 - Carmelo Mordini <carmelo> <cmordini@phys.ethz.ch>

import numpy as np
from pytrans.electrode import DCElectrode, RFElectrode
from pytrans.abstract_model import AbstractTrapModel
from rectset import rectangle_electrode as rect, pseudopotential as ps
from .geometry import corners, dc_width, rf_sep, rf_width

import logging
logger = logging.getLogger(__name__)

_DC_ELECTRODES = [
    "DCintop", "DCinbot",
    "DCtop1", "DCtop2", "DCtop3", "DCtop4", "DCtop5",
    "DCbot1", "DCbot2", "DCbot3", "DCbot4", "DCbot5",
]


class SurfaceTrapDCElectrode(DCElectrode):

    def __init__(self, name):
        self._corners = corners[name]
        super().__init__()

    def _unit_potential(self, x, y, z):
        return rect.rect_el_potential(x, y, z, *self._corners)

    def _unit_gradient(self, x, y, z):
        return rect.rect_el_gradient(x, y, z, *self._corners)

    def _unit_hessian(self, x, y, z):
        return rect.rect_el_hessian(x, y, z, *self._corners)


class SurfaceTrapRFElectrode(RFElectrode):

    def _unit_potential(self, x, y, z):
        return ps.pseudo_potential(x, y, z, rf_sep, rf_width)

    def _unit_gradient(self, x, y, z):
        return ps.pseudo_gradient(x, y, z, rf_sep, rf_width)

    def _unit_hessian(self, x, y, z):
        return ps.pseudo_hessian(x, y, z, rf_sep, rf_width)


class SurfaceTrap(AbstractTrapModel):

    _dc_electrodes = {name: SurfaceTrapDCElectrode(name) for name in _DC_ELECTRODES}
    _rf_electrodes = {"RF": SurfaceTrapRFElectrode()}
    _rf_voltage = 40
    _rf_freq_mhz = 20

    # Extra attributes and methods to enrich the model
    w_ele = dc_width
    x = np.arange(-350, 351, 0.5) * 1e-6
    y0 = 0.0
    z0 = ps.rf_null_z(rf_sep, rf_width)
    rf_null_coords = (None, y0, z0)
    dt = 392e-9

    @classmethod
    def x_ele(cls, j):
        # center position of electrode j
        # j = 1 .. 2
        return cls.w_ele * (j - 3)
