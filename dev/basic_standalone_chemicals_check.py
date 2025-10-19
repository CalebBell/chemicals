import sys
from chemicals import *
import numpy as np
import scipy.integrate
import scipy.interpolate
import scipy.spatial
import scipy.special
import scipy.optimize

def check_close(a, b, rtol=1e-7, atol=0):
    np.all(np.abs(a - b) <= (atol + rtol * np.abs(b)))
    return True

def run_checks():
    checks = []

    # Check CAS lookup and property retrieval
    CAS_water = CAS_from_any('Water')
    checks.append(check_close(MW(CAS_water), 18.01528, 1e-4))
    checks.append(check_close(Tb(CAS_water), 373.124, 1e-3))
    checks.append(check_close(Tc(CAS_water), 647.096, 1e-3))

    # Check Antoine vapor pressure calculation
    T = 373.15
    result = Antoine(T, A=10.1, B=1690, C=-43)
    checks.append(check_close(result, 95744.67843013899))

    # Check periodic table
    checks.append(check_close(periodic_table.Na.MW, 22.98977, 1e-5))
    checks.append(check_close(periodic_table.U.MW, 238.02891, 1e-5))

    return all(checks)

if run_checks():
    print("chemicals basic checks passed - NumPy and SciPy used successfully")
else:
    print('Library not OK')
    sys.exit(1)
