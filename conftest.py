import sys
import platform
is_pypy = 'PyPy' in sys.version
ver_tup = platform.python_version_tuple()[0:2]

def pytest_ignore_collect(path):
    path = str(path)
    if 'benchmark' in path or 'manual_runner' in path or 'make_test_stubs' in path or 'plot' in path or 'prerelease' in path:
        return True
    if 'dev' in path:
        return True 
#    if 'numba' in path:
#        return True
#    if 'units' in path:
#        return True
    if 'setup' in path:
        return True
    if ver_tup < ('3', '6') or ver_tup >= ('3', '9') or is_pypy:
        # numba does not yet run under pypy, numba support 3.6 to 3.8 for now
        if 'numba' in path:
            return True
