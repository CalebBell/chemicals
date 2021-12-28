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
from math import isnan, nan
try:
    path_join = os.path.join
except: # pragma: no cover
    pass
try:
    import numpy as np
    object_dtype = np.dtype(object)
    float64_dtype = np.dtype(np.float64)
    float32_dtype = np.dtype(np.float32)
    int64_dtype = np.dtype(np.int64)
    int32_dtype = np.dtype(np.int32)
    int16_dtype = np.dtype(np.int16)
    int8_dtype = np.dtype(np.int8)
    float_dtype_ids = set([id(float64_dtype), id(float32_dtype)])
    int_dtype_ids = set([id(int64_dtype), id(int32_dtype), id(int16_dtype), id(int8_dtype)])
    
except:
    pass
from chemicals.identifiers import CAS_to_int
# %% Loading data from local databanks

pd = None

df_sources = {}
load_cmds = {}

def make_df_sparse(df, non_sparse_columns=[]):
    '''Take a dataframe, and convert any floating-point columns which are mostly
    missing into sparse series. Return the resulting dataframe.
    '''
    sparse_float = pd.SparseDtype("float", nan)
    for col, dtype in zip(df.columns, df.dtypes):
        if col not in non_sparse_columns and id(dtype) in float_dtype_ids:
            series_orig = df[col]
            series_small = series_orig.astype(sparse_float)
            if series_small.memory_usage() < series_orig.memory_usage():
                df[col] = series_small
    return df


def register_df_source(folder, name, sep='\t', index_col=0, csv_kwargs={},
                       postload=None, sparsify=False, int_CAS=False):
    load_cmds[name] = (folder, name, sep, index_col, csv_kwargs, postload, sparsify, int_CAS)

'''The following flags will strip out the excess memory usage of redundant 
chemical metadata information.
'''
try:
    low_mem = bool(int(os.environ.get('CHEDL_LOW_MEMORY', '0')))
except:
    low_mem = False

spurious_columns = set(['name', 'formula', 'MW', 'InChI', 'InChI_key', 'Chemical',
                    'Data Type', 'Uncertainty', 'Fluid', 'Name', 'Names', 'Name ',
                    'Formula', 'Formula '])

def load_df(key):
    global pd
    if pd is None:
        import pandas as pd
    folder, name, sep, index_col, csv_kwargs, postload, sparsify, int_CAS = load_cmds[key]
    path = path_join(folder, name)
    df = pd.read_csv(path, sep=sep, index_col=index_col, **csv_kwargs)
    if postload: postload(df)
    if sparsify:
        df = make_df_sparse(df)
    if low_mem:
        for col_name in df.columns.values.tolist():
            if col_name in spurious_columns:
                del df[col_name]
                
    if int_CAS:
        df.index = pd.Index([CAS_to_int(s) for s in df.index], dtype=int64_dtype)
        
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
    except TypeError: # pragma: no cover
        raise TypeError("Method must be a string, not a %s object" %(type(method).__name__))
    return retrieve_from_df(df, index, key)

def retrieve_any_from_df_dict(df_dict, index, key):
    for df in df_dict.values():
        value = retrieve_from_df(df, index, key)
        if value is not None: return value

def retrieve_from_df(df, index, key):
    df_index = df.index
    if df_index.dtype is int64_dtype and isinstance(index, str):
        index = CAS_to_int(index)
    if index in df_index:
        if isinstance(key, (int, str)):
            return get_value_from_df(df, index, key)
        else: # Assume its an iterable of strings
            return [float(df.at[index, i]) for i in key]

def retrieve_any_from_df(df, index, keys):
    df_index = df.index
    if df_index.dtype is int64_dtype and isinstance(index, str):
        index = CAS_to_int(index)
    if index not in df.index: return None
    for key in keys:
        value = df.at[index, key]
        if not isnan(value):
            try:
                return float(value)
            except: # pragma: no cover
                return value

def get_value_from_df(df, index, key):
    value = df.at[index, key]
    try:
        return None if isnan(value) else float(value)
    except TypeError:
        # Not a number
        return value

def list_available_methods_from_df_dict(df_dict, index, key):
    methods = []
    for method, df in df_dict.items():
        df_index = df.index
        if df_index.dtype is int64_dtype and isinstance(index, str):
            index = CAS_to_int(index)
        if (index in df_index) and not isnan(df.at[index, key]):
            methods.append(method)
    return methods

def list_available_methods_from_df(df, index, keys_by_method):
    if index in df.index:
        return [method for method, key in keys_by_method.items()
                if not pd.isnull(df.at[index, key])]
    else:
        return []