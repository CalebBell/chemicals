"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2021, Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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
SOFTWARE.
"""

import numpy as np
import pytest

import chemicals
from chemicals.identifiers import CAS_to_int, check_CAS, int_to_CAS


@pytest.mark.slow
def test_CAS_numbers_valid_and_unique():
    int64_dtype = np.dtype(np.int64)
    chemicals.complete_lazy_loading()
    for k, df in chemicals.data_reader.df_sources.items():
        if df.index.dtype is int64_dtype:
            already_int = True
        else:
            already_int = False
        assert df.index.is_unique

        if already_int:
            for CAS in df.index.values.tolist():
                assert CAS < 9223372036854775807
                CAS = int_to_CAS(CAS)
                assert check_CAS(CAS)
        else:
            for CAS in df.index.values.tolist():
                assert check_CAS(CAS)
                CAS_int = CAS_to_int(CAS)
                # Check that the CAS number fits in a 64 bit int
                assert CAS_int < 9223372036854775807

    # Check that the name is CAS
    for k, df in chemicals.data_reader.df_sources.items():
        assert df.index.name == 'CAS'
