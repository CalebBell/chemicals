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


replaced = {'sum': np.sum}
replaced, NUMERICS_SUBMOD = fluids.numba.create_numerics(replaced, vec=False)


fluids.numba.transform_module(normal, __funcs, replaced, vec=False)


globals().update(__funcs)
globals().update(replaced)

def return_value_numpy(source):
    ret = re.search(r'return +\[', source)
    if ret:
        start_return, start_bracket = ret.regs[-1]
        enclosing = 1
        for i, v in enumerate(source[start_bracket:]):
            if v == '[':
                enclosing += 1
            if v == ']':
                enclosing -= 1
            if not enclosing:
                break
        return source[:start_bracket-1] + 'np.array([%s)' %source[start_bracket:i+start_bracket+1]        
    return source


# Magic to make a lists into arrays
list_mult_expr = r'\[ *([+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)) *\] *\* *([a-zA-Z0-9]+)'
numpy_not_list_expr = r'np.full((\4,), \1)'

import inspect
to_change = ['utils.zs_to_ws', 'utils.ws_to_zs', 'utils.zs_to_Vfs',
             'utils.dxs_to_dxsn1', 'utils.dxs_to_dns', 'utils.dns_to_dn_partials',
             'utils.dxs_to_dn_partials', 'utils.dxs_to_dxsn1',
             'utils.d2xs_to_dxdn_partials']
for s in to_change:
    mod, func = s.split('.')
    source = inspect.getsource(getattr(getattr(chemicals, mod), func))
    source = return_value_numpy(source)
#    source = source.replace(', kwargs={}', '').replace(', **kwargs', '')
    source = re.sub(list_mult_expr, numpy_not_list_expr, source)
    exec(source, globals(), globals())

    obj = numba.jit(cache=False)(globals()[func])
    globals()[func] = obj
    obj.__doc__ = ''
    
#    setattr(globals()[mod], func, obj)
    
    
    
    
    



'''
Add Function to use regular expressions, replace [0.0]*N by np.zeros
Will make functions like zs_to_ws work smoothly at optimal performance.


'''
