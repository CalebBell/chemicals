import os
import platform
import sys

is_pypy = 'PyPy' in sys.version
ver_tup = platform.python_version_tuple()[0:2]
ver_tup = tuple(int(i) for i in ver_tup)
is_x86_or_x86_64 = platform.machine().lower() in ('i386', 'i686', 'x86', 'x86_64', 'amd64')

def pytest_ignore_collect(path):
    path = str(path)
    if 'benchmark' in path or 'manual_runner' in path or 'make_test_stubs' in path or 'plot' in path or 'prerelease' in path or 'conf.py' in path:
        return True
    if 'dev' in path:
        return True
    if 'conf.py' in path:
        return True
#    if 'numba' in path:
#        return True
#    if 'units' in path:
#        return True
    if 'setup' in path:
        return True
    if ver_tup <= (3, 6) or ver_tup >= (3, 13) or is_pypy or not is_x86_or_x86_64:
        # numba does not yet run under pypy, numba support 3.6 to 3.9 for now
        if 'numba' in path:
            return True

def pytest_addoption(parser):
    parser.addoption(
        "--enable-low-memory", action="store", default="0", help="my option: 0 or 1"
    )

def pytest_configure(config):
    os.environ["CHEDL_LOW_MEMORY"] = config.getoption("--enable-low-memory")

