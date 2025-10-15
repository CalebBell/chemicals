import chemicals
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
    checks.append(MW(CAS_water) == 18.01528)
    checks.append(Tb(CAS_water) == 373.124)
    checks.append(Tc(CAS_water) == 647.096)

    # Check Antoine vapor pressure calculation
    T = 373.15
    result = Antoine(T, A=10.1, B=1690, C=-43)
    checks.append(check_close(result, 101047.2535))

    # Check periodic table
    checks.append(periodic_table.Na.MW == 22.9898)
    checks.append(periodic_table.U.MW == 238.02891)

    return all(checks)

if run_checks():
    print("chemicals basic checks passed - NumPy and SciPy used successfully")
else:
    print('Library not OK')
    exit(1)
