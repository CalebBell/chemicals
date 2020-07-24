import sys

def pytest_ignore_collect(path):
    path = str(path)
    if 'manual_runner' in path or 'make_test_stubs' in path:
        return True
