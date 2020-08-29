# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, Caleb Bell <Caleb.Andrew.Bell@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from fluids.numerics import assert_close, assert_close1d, linspace
import pytest
import pandas as pd
import numpy as np
from chemicals.permittivity import *
from chemicals.permittivity import permittivity_data_CRC

def test_permittivity_data():
    assert permittivity_data_CRC.index.is_unique
    assert permittivity_data_CRC.shape == (1303, 9)


def test_permittivity_IAPWS():
    Ts = [238., 256., 273., 273., 323., 323., 373., 373., 510., 523., 614., 647., 673., 673., 773., 773., 873.]
    rhos = [975.06, 995.25, 999.83, 1180., 988.10,  1258., 958.46, 1110., 15.832, 900., 94.29, 358., 100., 900., 100., 900., 450.]
    permittivity_calc = [permittivity_IAPWS(T, rho) for T, rho in zip(Ts, rhos)]
    permittivity_exp = [106.31159697963018, 95.19633650530638, 87.96431108641572, 107.06291112337524, 69.96455480833566, 97.7606839686273, 55.56584297721836, 67.73206302035597, 1.1224589212024803, 32.23229227177932, 1.7702660877662086, 6.194373838662447, 1.7541419715602131, 23.59653593827129, 1.6554135047590008, 20.160449425540186, 6.283091796558804]
    assert_close1d(permittivity_calc, permittivity_exp)


