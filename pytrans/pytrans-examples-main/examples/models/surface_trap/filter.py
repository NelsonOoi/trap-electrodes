#!/usr/bin/env python
# -*- coding:utf-8 -*-
#
# Created: 06/2023
# Author: Carmelo Mordini <cmordini@phys.ethz.ch>

from .model import SurfaceTrap
from scipy import signal as sg

dt = SurfaceTrap.dt
fs = 1 / dt  # waveform sample rate

# two-pole RC filter
R1 = 200
C1 = 10e-9
R2 = 180
C2 = 33e-9

tau1 = R1 * C1
tau2 = R2 * C2

zpk_analog = [], [-1 / tau1, -1 / tau2], 1 / tau1 / tau2
zpk_dig = sg.bilinear_zpk(*zpk_analog, fs)
# b, a = sg.zpk2tf(*zpk_dig)

trap_filter = sg.dlti(*zpk_dig, dt=dt)
