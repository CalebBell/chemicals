# -*- coding: utf-8 -*-
'''Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2020, Caleb Bell <Caleb.Andrew.Bell@gmail.com>
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
SOFTWARE.'''

from __future__ import division
import sys
import importlib.util
import types
import numpy as np
import numba
import chemicals
import fluids.numba
import re
normal = chemicals
import inspect

'''
Known to not work - reason
elements.py - everything is dict-based; separate function implementations?
viscosity.viscosity_converter - scipy splines
lennard_jones.collision_integral_Kim_Monroe - dict lookup
lennard_jones.collision_integral_Neufeld_Janzen_Aziz - dict lookup
temperature.T_converter - global, function-in-function, scipy splines
critical.Ihmels and friends - assert, adding None (consider re-designing interface, will make pint easier)
'''
__all__ = []

__funcs = {}


replaced = {'sum': np.sum, 'combinations': fluids.numba.combinations}
replaced, NUMERICS_SUBMOD = fluids.numba.create_numerics(replaced, vec=False)


blacklist = set(['to_num'])
fluids.numba.transform_module(normal, __funcs, replaced, vec=False,
                              blacklist=blacklist)



to_change = ['utils.zs_to_ws', 'utils.ws_to_zs', 'utils.zs_to_Vfs',
             'utils.dxs_to_dxsn1', 'utils.dxs_to_dns', 'utils.dns_to_dn_partials',
             'utils.dxs_to_dn_partials', 'utils.dxs_to_dxsn1',
             'utils.d2xs_to_dxdn_partials', 'viscosity.Lorentz_Bray_Clarke',
             'viscosity.Herning_Zipperer', 'volume.COSTALD_mixture',
             'rachford_rice.Rachford_Rice_solution', 'rachford_rice.Rachford_Rice_solution_LN2',
             'rachford_rice.Rachford_Rice_solution_numpy',
             'rachford_rice.Rachford_Rice_polynomial', 
             'rachford_rice.Rachford_Rice_polynomial_3',
             'rachford_rice.Rachford_Rice_polynomial_4',
             'rachford_rice.Rachford_Rice_polynomial_5',
             'rachford_rice.Rachford_Rice_solution_polynomial',
             'rachford_rice.Rachford_Rice_numpy_err_fprime2',
             'rachford_rice.Li_Johns_Ahmadi_solution',
             'rachford_rice._Rachford_Rice_analytical_3',
             'rachford_rice.flash_inner_loop',
             'rachford_rice.Rachford_Rice_solution2',
             'rachford_rice.Rachford_Rice_solutionN',
             'rachford_rice.RRN_new_betas',
             'rachford_rice.Rachford_Rice_flashN_f_jac',
             'rachford_rice.Rachford_Rice_flash2_f_jac',
             'rachford_rice.Rachford_Rice_valid_solution_naive',
             ]

fluids.numba.transform_lists_to_arrays(chemicals, to_change, __funcs)
globals().update(__funcs)
globals().update(replaced)