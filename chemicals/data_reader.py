# -*- coding: utf-8 -*-
'''Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2020 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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

__all__ = ['df_sources', 'data_source', 'register_df_source', 'load_df']

import os
import pandas as pd
path_join = os.path.join


df_sources = {}

load_cmds = {}

def register_df_source(folder, name, sep='\t', index_col=0, csv_kwargs={},
                       postload=None):
    load_cmds[name] = (folder, name, sep, index_col, csv_kwargs, postload)

def load_df(key):
    folder, name, sep, index_col, csv_kwargs, postload = load_cmds[key]
    path = path_join(folder, name)
    df = pd.read_csv(path, sep=sep, index_col=index_col, **csv_kwargs)
    if postload is not None:
        postload(df)
        
    df_sources[key] = df
    

def data_source(key):
    try:
        return df_sources[key]
    except KeyError:
        load_df(key)
        return df_sources[key]
