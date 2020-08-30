#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
try:
    import test_acentric
except:
    print('run this from the tests directory')
    exit()
import test_combustion
import test_critical
import test_dipole
import test_dippr
import test_elements
import test_environment
import test_heat_capacity
import test_interface
import test_lennard_jones
import test_miscdata
import test_permittivity
import test_phase_change
import test_rachford_rice
import test_reaction
import test_refractivity
import test_solubility
import test_temperature
import test_thermal_conductivity
import test_triple
import test_utils
import test_vapor_pressure
import test_vectorized
import test_virial
import test_viscosity
import test_volume

# dynamically generated code - numba, units, vectorize - not part of this test suite
to_test = [test_acentric, test_combustion, test_critical, test_dipole, test_dippr, test_elements, test_environment, test_heat_capacity, test_interface, test_lennard_jones, test_miscdata, test_permittivity, test_phase_change, test_rachford_rice, test_reaction, test_refractivity, test_solubility, test_temperature, test_thermal_conductivity, test_triple, test_utils, test_vapor_pressure, test_vectorized, test_virial, test_viscosity, test_volume]


skip_marks = ['slow', 'fuzz', 'skip_types']
skip_marks_set = set(skip_marks)
if len(sys.argv) >= 2:
    #print(sys.argv)
    # Run modules specified by user
    to_test = [globals()[i] for i in sys.argv[1:]]
for mod in to_test:
    print(mod)
    for s in dir(mod):
        skip = False
        obj = getattr(mod, s)
        if callable(obj) and hasattr(obj, '__name__') and obj.__name__.startswith('test'):
            try:
                for bad in skip_marks:
                    if bad in obj.__dict__:
                        skip = True
                if 'pytestmark' in obj.__dict__:
                    marked_names = [i.name for i in obj.__dict__['pytestmark']]
                    for mark_name in marked_names:
                        if mark_name in skip_marks_set:
                            skip = True
            except Exception as e:
                #print(e)
                pass
            if not skip:
                try:
                    print(obj)
                    obj()
                except Exception as e:
                    print('FAILED TEST %s with error:' %s)
                    print(e)
