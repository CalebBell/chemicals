# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL).

Utilities for process modeling. Copyright (C) 2020, Caleb Bell
<Caleb.Andrew.Bell@gmail.com> Permission is hereby granted, free of charge, to
any person obtaining a copy of this software and associated documentation files
(the "Software"), to deal in the Software without restriction, including without
limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom
the Software is furnished to do so, subject to the following conditions: The
above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS
IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from chemicals.utils import PY37, numba_blacklisted, numba_cache_blacklisted
import fluids as normal_fluids

busy = False
__all__ = []

to_change = ['utils.zs_to_ws', 'utils.ws_to_zs', 'utils.zs_to_Vfs',
     'utils.dxs_to_dxsn1', 'utils.dxs_to_dns', 'utils.dns_to_dn_partials',
     'utils.dxs_to_dn_partials', 'utils.dxs_to_dxsn1',
     'utils.d2xs_to_dxdn_partials', 'viscosity.Lorentz_Bray_Clarke',
     'viscosity.Herning_Zipperer', 'volume.COSTALD_mixture',
     'viscosity.Wilke_large', 'viscosity.Wilke_prefactored',
     'interface.Winterfeld_Scriven_Davis',

     'rachford_rice.Rachford_Rice_solution',
     'rachford_rice.Rachford_Rice_solution_LN2',
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
     'flash_basic.flash_wilson',
     'solubility.Henry_pressure_mixture',
     'critical.Chueh_Prausnitz_Tc',
     'critical.Chueh_Prausnitz_Vc',
     'thermal_conductivity.Lindsay_Bromley',
     'thermal_conductivity.DIPPR9I',
     'virial.Z_from_virial_density_form',
     
     'vapor_pressure.Wagner_fitting_jacobian',
     'vapor_pressure.Wagner_original_fitting_jacobian',
     'vapor_pressure.Antoine_fitting_jacobian',
     
     'dippr.EQ102_fitting_jacobian',
     'dippr.EQ101_fitting_jacobian',
     'dippr.EQ106_fitting_jacobian',
     'dippr.EQ105_fitting_jacobian',
     'dippr.EQ107_fitting_jacobian',
     'dippr.EQ102',
     
     'vapor_pressure.Yaws_Psat_fitting_jacobian',
     'viscosity.mu_Yaws_fitting_jacobian',
     ]

def transform_complete_chemicals(replaced, __funcs, __all__, normal, vec):
    cache_blacklist = set(numba_cache_blacklisted)
    __funcs.update(normal_fluids.numba.numbafied_fluids_functions.copy())
    new_mods = normal_fluids.numba.transform_module(normal, __funcs, replaced, vec=vec,
                                                    blacklist=set(numba_blacklisted),
                                                    cache_blacklist=cache_blacklist)

    normal_fluids.numba.transform_lists_to_arrays(normal, to_change, __funcs, cache_blacklist=cache_blacklist, vec=vec)

    for mod in new_mods:
        mod.__dict__.update(__funcs)
        try:
            __all__.extend(mod.__all__)
        except AttributeError:
            pass

def transform():
    gdct = globals()
    import chemicals
    import fluids.numba
    normal = chemicals
    __funcs = {}
    replaced = fluids.numba.numerics_dict.copy()
    orig_file = __file__
    transform_complete_chemicals(replaced, __funcs, __all__, normal, vec=False)
    gdct.update(__funcs)
    gdct.update(replaced)
    gdct['__name__'] = 'chemicals.numba'
    gdct['__file__'] = orig_file
    del gdct['transform']

if PY37:
    def __getattr__(name):
        if name == '__path__' or busy:
            raise AttributeError("module %s has no attribute %s" %(__name__, name))
        gdct = globals()
        gdct['busy'] = True
        transform()
        del gdct['__getattr__'], gdct['busy']
        try: return gdct[name]
        except KeyError:
            raise AttributeError("module %s has no attribute %s" %(__name__, name))
else:
    transform()
