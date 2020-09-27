import sys

def pytest_ignore_collect(path):
    path = str(path)
    if 'manual_runner' in path or 'make_test_stubs' in path:
        return True
    if 'dev' in path:
        return True 
    if 'numba' in path:
        return True
#    if 'units' in path:
#        return True
    if 'setup' in path:
        return True
