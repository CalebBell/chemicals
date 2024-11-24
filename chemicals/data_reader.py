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
    float_dtype_ids = {id(float64_dtype), id(float32_dtype)}
    int_dtype_ids = {id(int64_dtype), id(int32_dtype), id(int16_dtype), id(int8_dtype)}

except:
    pass
try:
    import threading
    import sqlite3
except:
    pass
from chemicals.identifiers import CAS_to_int
from chemicals.utils import source_path

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


def register_df_source(folder, name, sep='\t', index_col=0, csv_kwargs=None,
                       postload=None, sparsify=False, int_CAS=False):
    if csv_kwargs is None: csv_kwargs = {}
    load_cmds[name] = (folder, name, sep, index_col, csv_kwargs, postload, sparsify, int_CAS)

"""The following flags will strip out the excess memory usage of redundant
chemical metadata information.
"""
try:
    low_mem = bool(int(os.environ.get('CHEDL_LOW_MEMORY', '0')))
except:
    low_mem = False

spurious_columns = {'name', 'formula', 'MW', 'InChI', 'InChI_key', 'Chemical',
                    'Data Type', 'Uncertainty', 'Fluid', 'Name', 'Names', 'Name ',
                    'Formula', 'Formula '}

def load_df(key):
    global pd
    if pd is None:
        import pandas as pd
    folder, name, sep, index_col, csv_kwargs, postload, sparsify, int_CAS = load_cmds[key]
    path = path_join(folder, name)
    if int_CAS:
        dtype = csv_kwargs.get('dtype', {})
        dtype['CAS'] = int64_dtype
        csv_kwargs['dtype'] = dtype
    df = pd.read_csv(path, sep=sep, index_col=index_col, **csv_kwargs)
    if postload: postload(df)
    if sparsify:
        df = make_df_sparse(df)
    if low_mem:
        for col_name in df.columns.values.tolist():
            if col_name in spurious_columns:
                df[col_name] = pd.Series([], dtype=float).astype(pd.SparseDtype("float", nan))

    if int_CAS and df.index.dtype is object_dtype:
        # If the index is already an int, leave it be
        """Most CAS numbers fit in 32 bits. Not all of them do though, for
        example https://commonchemistry.cas.org/detail?cas_rn=2222298-66-8
        or 2627558-64-7

        The maximum value of an unsigned 32 bit integer is 4294967295.
        For 64 bit it is 18446744073709551615.

        It would be possible to remove the check digit of the CAS number,
        which would allow all 10-digit current CAS format integers to fit
        into an unsigned 32 bit integer.
        https://www.cas.org/support/documentation/chemical-substances/faqs
        CAS says they are only "up to ten digits". However, before 2008 all
        CAS numbers were "up to nine digits"; and they are already 25% of the
        way through 10 digits. It is only a matter of time before they add
        another digit. At their current rate this will be in 2036, but it will
        likely be well before that. Therefore, it does not justify removing
        the check digit.
        """
        df.index = pd.Index([CAS_to_int(s) for s in df.index], dtype=int64_dtype, name=df.index.name)

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
        raise ValueError(f'Invalid method: {method}, allowed methods are {list(df_dict)}')
    except TypeError: # pragma: no cover
        raise TypeError(f"Method must be a string, not a {type(method).__name__} object")
    return retrieve_from_df(df, index, key)

def retrieve_any_from_df_dict(df_dict, index, key):
    for df in df_dict.values():
        value = retrieve_from_df(df, index, key)
        if value is not None: return value

def retrieve_from_df(df, index, key):
    df_index = df.index
    if df_index.dtype is int64_dtype and isinstance(index, str):
        try: index = CAS_to_int(index)
        except: return None
    if index in df_index:
        if isinstance(key, (int, str)):
            return get_value_from_df(df, index, key)
        else: # Assume its an iterable of strings
            return [float(df.at[index, i]) for i in key]

def retrieve_any_from_df(df, index, keys):
    df_index = df.index
    if df_index.dtype is int64_dtype and isinstance(index, str):
        try: index = CAS_to_int(index)
        except: return None
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
    except TypeError: # Not a number
        return value

def list_available_methods_from_df_dict(df_dict, index, key):
    methods = []
    int_index = None if type(index) is str else index # Assume must be string or int
    for method, df in df_dict.items():
        df_index = df.index
        if df_index.dtype is int64_dtype:
            if int_index is None:
                int_index = CAS_to_int(index)
            if int_index in df_index and not isnan(df.at[int_index, key]):
                methods.append(method)
        elif (index in df_index) and not isnan(df.at[index, key]):
            methods.append(method)
    return methods

def list_available_methods_from_df(df, index, keys_by_method):
    if index in df.index:
        return [method for method, key in keys_by_method.items()
                if not pd.isnull(df.at[index, key])]
    else:
        return []

### Database

try:
    USE_CONSTANTS_DATABASE = os.path.exists(path_join(source_path, 'Misc', 'default.sqlite'))
except:
    USE_CONSTANTS_DATABASE = False

CONSTANT_DATABASE_COLUMNS = ['index', 'MW', 'Tt', 'Tm', 'Tb', 'Tc', 'Pt', 'Pc', 'Vc',
'Zc', 'omega', 'T_flash', 'T_autoignition', 'LFL', 'UFL',
'Hfs', 'Hfl', 'Hfg', 'S0s', 'S0l', 'S0g',
'Hfus', 'Stockmayer', 'molecular_diameter',
'dipole_moment', 'logP', 'RG', 'RON', 'MON', 'ignition_delay', 'linear',
'GWP', 'ODP', 'RI', 'RIT', 'GTP', 'HANSEN_DELTA_D', 'HANSEN_DELTA_P', 'HANSEN_DELTA_H']



CONSTANT_DATABASE_NAME_TO_IDX = {k: i for i, k in enumerate(CONSTANT_DATABASE_COLUMNS)}
CONSTANTS_CURSOR = None
try:
    thread_local_storage = threading.local()
except:
    pass

DATABASE_CONSTANTS_CACHE = {}

def cached_constant_lookup(CASi, prop):
    '''Look up a constant property for a compound, either from cache or the database.'''
    if not hasattr(thread_local_storage, 'cursor'):
        init_constants_db()

    # Check if the result is cached
    if CASi in DATABASE_CONSTANTS_CACHE:
        result = DATABASE_CONSTANTS_CACHE[CASi]
    else:
        # Fetch the whole row from the database
        thread_local_storage.cursor.execute("SELECT * FROM constants WHERE `index`=?", (str(CASi),))
        result = thread_local_storage.cursor.fetchone()
        DATABASE_CONSTANTS_CACHE[CASi] = result

    if result is None:
        return result, False  # Return result and indicate the compound was not found

    # Retrieve the value for the specified property
    prop_idx = CONSTANT_DATABASE_NAME_TO_IDX[prop]
    return result[prop_idx], True

def init_constants_db():
    '''Initialize the database connection and cursor for the current thread if not already done.'''
    if not hasattr(thread_local_storage, 'conn'):
        # Create a new connection and cursor for the thread
        thread_local_storage.conn = sqlite3.connect(
            path_join(source_path, 'Misc', 'default.sqlite'),
            check_same_thread=False
        )
        thread_local_storage.cursor = thread_local_storage.conn.cursor()

def database_constant_lookup(CASi, prop):
    if type(CASi) is str: # Assume it must be either an int or string
        try:
            CASi = CAS_to_int(CASi)
        except:
            return None, False
    try:
        return cached_constant_lookup(CASi, prop)
    except (TypeError, KeyError) as e:
        raise e
    except Exception:
        # Prevent database lookup after first failure considering it should work everytime.
        # It will possibly fail for users every time if database has not been created.
        global USE_CONSTANTS_DATABASE
        USE_CONSTANTS_DATABASE = False
        return None, False
