# -*- coding: utf-8 -*-
"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
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
SOFTWARE.
"""
__all__ = ['df_sources', 
           'data_source', 
           'register_df_source', 
           'load_df',
           'retrieve_any_from_df_dict',
           'retrieve_from_df_dict',
           'retrieve_any_from_df',
           'retrieve_from_df',
           'list_available_methods_from_df_dict',
           'list_available_methods_from_df']

import os
from math import isnan
try:
    from collections.abc import Iterable
except:
    try:
        from collections import Iterable
    except:
        Iterable = list
try:
    path_join = os.path.join
except:
    pass

# %% Loading data from local databanks

pd = None

df_sources = {}
load_cmds = {}

def register_df_source(folder, name, sep='\t', index_col=0, csv_kwargs={},
                       postload=None):
    load_cmds[name] = (folder, name, sep, index_col, csv_kwargs, postload)

def load_df(key):
    global pd
    if pd is None:
        import pandas as pd
    folder, name, sep, index_col, csv_kwargs, postload = load_cmds[key]
    path = path_join(folder, name)
    df = pd.read_csv(path, sep=sep, index_col=index_col, **csv_kwargs)
    if postload: postload(df)
    df_sources[key] = df

def data_source(key):
    try:
        return df_sources[key]
    except KeyError:
        load_df(key)
        return df_sources[key]


# %% Retrieving data from files

def retrieve_from_df_dict(df_dict, index, key, method):
    try: 
        df = df_dict[method]
    except KeyError: 
        raise ValueError('Invalid method: %s, allowed methods are %s' %(
                method, list(df_dict)))
    except TypeError: 
        raise TypeError("Method must be a string, not a %s object" %(type(method).__name__))
    return retrieve_from_df(df, index, key)

def retrieve_any_from_df_dict(df_dict, index, key):
    for df in df_dict.values():
        value = retrieve_from_df(df, index, key)
        if value is not None: return value

def retrieve_from_df(df, index, key):
    if index in df.index:
        if isinstance(key, str):
            return get_value_from_df(df, index, key)
        elif isinstance(key, Iterable):
            return [float(df.at[index, i]) for i in key]
        else:
            raise ValueError('key must be a string or an iterable of strings')

def retrieve_any_from_df(df, index, keys):
    if isinstance(keys, str) or not isinstance(keys, Iterable):    
        raise ValueError('keys must be an iterable of strings')
    if index not in df.index: return None
    for key in keys:
        value = df.at[index, key]
        if not isnan(value):
            try:
                return float(value)
            except:
                return value

def get_value_from_df(df, index, key):
    value = df.at[index, key]
    try:
        return None if isnan(value) else float(value)
    except TypeError:
        # Not a number
        return value
            
def list_available_methods_from_df_dict(df_dict, index, key):
    return [method for method, df in df_dict.items()
            if (index in df.index) and not pd.isnull(df.at[index, key])]

def list_available_methods_from_df(df, index, keys):
    if index in df.index:
        return [key for key in keys if not pd.isnull(df.at[index, key])]
    else:
        return []