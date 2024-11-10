"""Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019, 2020, 2021, 2022 Caleb Bell
<Caleb.Andrew.Bell@gmail.com>

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


import fluids

__version__ = '1.3.1'
from math import isnan

if not fluids.numerics.is_micropython:

    from . import (
        acentric,
        air,
        combustion,
        critical,
        data_reader,
        dipole,
        dippr,
        elements,
        environment,
        exceptions,
        flash_basic,
        heat_capacity,
        iapws,
        identifiers,
        interface,
        lennard_jones,
        miscdata,
        molecular_geometry,
        permittivity,
        phase_change,
        rachford_rice,
        reaction,
        refractivity,
        safety,
        solubility,
        temperature,
        thermal_conductivity,
        triple,
        utils,
        vapor_pressure,
        virial,
        viscosity,
        volume,
    )
    from .acentric import *  # noqa: F403
    from .air import *  # noqa: F403
    from .combustion import *  # noqa: F403
    from .critical import *  # noqa: F403
    from .dipole import *  # noqa: F403
    from .dippr import *  # noqa: F403
    from .elements import *  # noqa: F403
    from .environment import *  # noqa: F403
    from .exceptions import *  # noqa: F403
    from .flash_basic import *  # noqa: F403
    from .heat_capacity import *  # noqa: F403
    from .iapws import *  # noqa: F403
    from .identifiers import *  # noqa: F403
    from .interface import *  # noqa: F403
    from .lennard_jones import *  # noqa: F403
    from .miscdata import *  # noqa: F403
    from .molecular_geometry import *  # noqa: F403
    from .permittivity import *  # noqa: F403
    from .phase_change import *  # noqa: F403
    from .rachford_rice import *  # noqa: F403
    from .reaction import *  # noqa: F403
    from .refractivity import *  # noqa: F403
    from .safety import *  # noqa: F403
    from .solubility import *  # noqa: F403
    from .temperature import *  # noqa: F403
    from .thermal_conductivity import *  # noqa: F403
    from .triple import *  # noqa: F403
    from .utils import *  # noqa: F403
    from .utils import PY37
    from .vapor_pressure import *  # noqa: F403
    from .virial import *  # noqa: F403
    from .viscosity import *  # noqa: F403
    from .volume import *  # noqa: F403
    __all__ = ['utils', 'critical', 'elements', 'reaction', 'dipole', 'dippr',
               'temperature', 'miscdata', 'environment', 'refractivity', 'solubility',
               'lennard_jones', 'heat_capacity', 'vapor_pressure', 'virial',
               'phase_change', 'triple', 'exceptions', 'acentric', 'viscosity',
               'interface', 'permittivity', 'thermal_conductivity', 'combustion',
               'volume', 'rachford_rice', 'flash_basic', 'identifiers', 'safety',
               'iapws', 'air', 'molecular_geometry']

    __all__.extend(exceptions.__all__)
    __all__.extend(critical.__all__)
    __all__.extend(utils.__all__)
    __all__.extend(elements.__all__)
    __all__.extend(reaction.__all__)
    __all__.extend(dipole.__all__)
    __all__.extend(dippr.__all__)
    __all__.extend(temperature.__all__)
    __all__.extend(miscdata.__all__)
    __all__.extend(environment.__all__)
    __all__.extend(refractivity.__all__)
    __all__.extend(solubility.__all__)
    __all__.extend(lennard_jones.__all__)
    __all__.extend(heat_capacity.__all__)
    __all__.extend(vapor_pressure.__all__)
    __all__.extend(virial.__all__)
    __all__.extend(phase_change.__all__)
    __all__.extend(triple.__all__)
    __all__.extend(acentric.__all__)
    __all__.extend(viscosity.__all__)
    __all__.extend(interface.__all__)
    __all__.extend(permittivity.__all__)
    __all__.extend(thermal_conductivity.__all__)
    __all__.extend(combustion.__all__)
    __all__.extend(volume.__all__)
    __all__.extend(rachford_rice.__all__)
    __all__.extend(flash_basic.__all__)
    __all__.extend(identifiers.__all__)
    __all__.extend(safety.__all__)
    __all__.extend(iapws.__all__)
    __all__.extend(air.__all__)
    __all__.extend(molecular_geometry.__all__)

    submodules = [critical, utils, elements, dipole, dippr, temperature, miscdata,
                  environment, refractivity, solubility, lennard_jones,
                  vapor_pressure, virial, phase_change, triple, acentric, viscosity,
                  interface, permittivity, thermal_conductivity, combustion,
                  heat_capacity, reaction, volume, rachford_rice, flash_basic, identifiers, safety,
                  iapws, air, molecular_geometry]


    chemicals_dir = utils.source_path


    def complete_lazy_loading():
        critical._load_critical_data()
        dipole._load_dipole_data()
        environment._load_GWP_ODP_data()
        environment._load_logP_data()
        heat_capacity._load_Cp_data()
        interface.load_interface_dfs()
        lennard_jones._load_LJ_data()
        miscdata._load_VDI_saturation_dict()
        miscdata._load_miscdata()
        permittivity._load_permittivity_data()
        phase_change._load_phase_change_constants()
        phase_change._load_phase_change_correlations()
        reaction._load_reaction_data()
        refractivity._load_RI_data()
        safety._load_safety_data()
        thermal_conductivity._load_k_data()
        triple._load_triple_data()
        vapor_pressure.load_vapor_pressure_dfs()
        viscosity._load_mu_data()
        volume._load_rho_data()
        molecular_geometry._load_RG_data()
        combustion._load_combustion_data()
        solubility._load_solubility_data()
        try:
            identifiers.search_chemical('asdfasddsaf', autoload=True, cache=False)
        except:
            pass

    def remove_missing(values, CASs):
        values_found, CASs_found = [], []
        for i in range(len(values)):
            if not isnan(values[i]):
                values_found.append(values[i])
                CASs_found.append(CASs[i])
        return values_found, CASs_found
    def constant_statistics():
        '''

        Returns
        -------
        dict : str[int]
            Counts of each constant property and the number of chemicals that
            values are available for, [-]

        '''
        complete_lazy_loading()
        counts = {}
        properties = ['Tt', 'Pt', 'Tm', 'Tb', 'Tc', 'Pc', 'Vc', 'Zc',
                      'omega', 'T_flash', 'T_autoignition', 'LFL', 'UFL',
                     'Hfs', 'Hfl', 'Hfg', 'S0s', 'S0l', 'S0g',
                     'RI', 'Hfus', 'Stockmayer', 'molecular_diameter'
                     'Dipole', 'logP', 'RG', 'RON', 'MON', 'IGNITION_DELAY',
                     'linear', 'HANSEN_DELTA_D', 'HANSEN_DELTA_P', 'HANSEN_DELTA_H']
        for df in data_reader.df_sources.values():
            for p in properties:
                if p in df.columns:
                    prop_values, prop_CASs = df[p].values, df.index.values
                    prop_values, prop_CASs = remove_missing(prop_values, prop_CASs)
                    if p in counts:
                        counts[p].update(prop_CASs)
                    else:
                        counts[p] = set(prop_CASs)
        return {k: len(v) for k, v in counts.items()}

    def memory_usage(finish_loading=True):
        if finish_loading:
            complete_lazy_loading()
        usages = []
        names = []
        for name, df in data_reader.df_sources.items():
            names.append(name)
            usages.append(df.memory_usage().values.sum())

        names = [x for _, x in sorted(zip(usages, names))]
        usages.sort()
        for name, use in zip(names, usages):
            print(f'{name} : {use/1024**2:3f} MB')
        print(f'Total usage: {sum(usages)/1024**2:3f} MB')

    if PY37:
        def __getattr__(name):
            global vectorized, numba, units, numba_vectorized
            if name == 'vectorized':
                import chemicals.vectorized
                return chemicals.vectorized
            if name == 'numba':
                import chemicals.numba
                return chemicals.numba
            if name == 'units':
                import chemicals.units
                return chemicals.units
            if name == 'numba_vectorized':
                import chemicals.numba_vectorized
                return chemicals.numba_vectorized
            raise AttributeError(f"module {__name__} has no attribute {name}")
    else:
        pass
