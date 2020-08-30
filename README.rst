=========
chemicals
=========

.. image:: http://img.shields.io/pypi/v/chemicals.svg?style=flat
   :target: https://pypi.python.org/pypi/chemicals
   :alt: Version_status
.. image:: http://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat
   :target: https://chemicals.readthedocs.io/en/latest/
   :alt: Documentation
.. image:: https://travis-ci.com/CalebBell/chemicals.svg?branch=release?style=flat
   :target: https://travis-ci.com/github/CalebBell/chemicals
   :alt: Build_status
.. image:: http://img.shields.io/badge/license-MIT-blue.svg?style=flat
   :target: https://github.com/CalebBell/chemicals/blob/master/LICENSE.txt
   :alt: license
.. image:: https://img.shields.io/coveralls/CalebBell/chemicals.svg
   :target: https://coveralls.io/github/CalebBell/chemicals
   :alt: Coverage
.. image:: https://img.shields.io/pypi/pyversions/chemicals.svg
   :target: https://pypi.python.org/pypi/chemicals
   :alt: Supported_versions
.. image:: https://badges.gitter.im/CalebBell/chemicals.svg
   :alt: Join the chat at https://gitter.im/CalebBell/chemicals
   :target: https://gitter.im/CalebBell/chemicals
.. image:: http://img.shields.io/appveyor/ci/calebbell/chemicals.svg
   :target: https://ci.appveyor.com/project/calebbell/chemicals/branch/master
   :alt: Build_status
.. image:: https://zenodo.org/badge/264697738.svg
   :alt: Zenodo
   :target: https://zenodo.org/badge/latestdoi/264697738

.. contents::

Key Features & Capabilities
---------------------------

The chemicals library features an extensive compilation of pure component 
chemical data that can serve engineers, scientists, technicians, and anyone 
working with chemicals. The chemicals library facilitates the retrieval of:

- Chemical constants including formula, molecular weight, normal boiling and 
  melting points, triple point, heat of formation, absolute entropy of 
  formation, heat of fusion, similarity variable, dipole moment, accentric 
  factor, etc.
  
- Assorted information of safety and toxicity of chemicals.

- Methods (and their respective coefficients) for the calculation of temperature
  and pressure dependent chemical properties including vapor pressure, 
  heat capacity, molar volume, thermal conductivity, surface tension, dynamic 
  viscosity, heat of vaporization, relative permitivity, etc.

Data for over 20,000 chemicals are made available as local databanks in this 
library. All databanks are loaded on-demand, saving loading time and RAM. For
example, if only data on the normal boiling point is required, the chemicals 
library will only load normal boiling point datasets. This on-demand loading 
feature makes the chemicals library an attractive dependence for software
modeling chemical processes. In fact, `The Biorefinery Simulation and Techno-Economic Analysis Modules (BioSTEAM) <https://biosteam.readthedocs.io/en/latest/>`_ 
is reliant on the chemicals library for the simulation of unit operations.

If you need to know something about a chemical, give chemicals a try.

# TODO: 

# Do we want to keep the title of this section?

# Should we prioritize a logo for ChEDL and have it in all libraries?

# Note how thermo is reliant on the chemicals library too

# Note on flash capabilities?

# Note on numba?

# Note on pint?

# Once we do finish working on the README, we should update the docs/index.rst file

Installation
------------

Get the latest version of chemicals from
https://pypi.python.org/pypi/chemicals/

If you have an installation of Python with pip, simple install it with:

    $ pip install chemicals
    
Alternatively, if you are using `conda <https://conda.io/en/latest/>`_ as your 
package management, you can install chemicals in your environment from 
`conda-forge <https://conda-forge.org/>`_ channel with:

    $ conda install -c conda-forge chemicals

To get the git version, run:

    $ git clone git://github.com/CalebBell/chemicals.git

# TODO: Add conda installation

Documentation
-------------

chemicals's documentation is available on the web:

    http://chemicals.readthedocs.io/

Getting Started
---------------

The library is designed around SI units. The retrieval of constant 
chemical properties is done by CASRN: 

.. code-block:: python

    >>> from chemicals import CAS_from_any, MW, Tb, Tm, Tc, Pc, Vc, Hfus, Hfs, Hfl, Hfg, S0s, S0l, S0g
    >>> # Search for CASRN by a common name or a unique identifier such as the IchI key
    >>> CAS_water = CAS_from_any('Water')
    >>> MW(CAS_water)   # Molecular weight [g/mol]
    18.01528
    >>> Tb(CAS_water)   # Normal boiling point [K]
    373.124
    >>> Tm(CAS_water)   # Melting point [K]
    273.15
    >>> Tc(CAS_water)   # Critical temperature [K]
    647.14
    >>> Pc(CAS_water)   # Critical pressure [Pa]
    22048320.0
    >>> Vc(CAS_water)   # Critical volume [m^3/mol]
    5.6000000000000006e-05
    >>> Hfus(CAS_water) # Heat of fusion [J/mol]
    6010.0
    >>> Hfs(CAS_water)  # Solid heat of formation [J/mol] 
    None 
    >>> Hfl(CAS_water)  # Liquid heat of formation [J/mol]
    -285825.0
    >>> Hfg(CAS_water)  # Gas heat of formation [J/mol] 
    -241822.0
    >>> S0s(CAS_water)  # Absolute solid enthalpy of formation [J/mol/K]
    None
    >>> S0l(CAS_water)  # Absolute liquid enthalpy of formation [J/mol/K]
    70.0
    >>> S0g(CAS_water)  # Absolute gas enthalpy of formation [J/mol/K]
    188.8

Such "lookup functions" search through many databanks until the value for the 
given CASRN is found. When no value is available, None is returned. Note that some 
properties, including the solid heat of formation (Hfs) and absolute entropy of 
formation (S0s) are missing for water. You can optionally pass a "method" to select 
which from databank to retrieve data:

.. code-block:: python

    >>> Tb(CAS_water, method='YAWS')
    373.15
 
To view all available methods for a given chemical, just use the functions with 
"_methods" tagged at the end of the name:

.. code-block:: python
    
    >>> Tb_methods(CAS_water)
    ['CRC_INORG', 'YAWS']

The databanks can also be accessed through their respective module:

>>> from chemicals.critical import critical_data_Yaws
>>> critical_data_Yaws
                               Chemical       Tc  ...  omega        Zc
CASRN                                             ...                 
100-00-5           p-chloronitrobenzene   751.00  ...  0.491  0.275355
100-01-6                 p-nitroaniline   851.00  ...  0.782  0.253621
100-10-7    p-dimethylaminobenzaldehyde   832.00  ...  0.527  0.209027
100-18-5           p-diisopropylbenzene   689.00  ...  0.390  0.255749
100-21-0              terephthalic acid  1113.00  ...  1.059  0.180981
                                ...      ...  ...    ...       ...
99814-65-0            1-pentadeceN-3-ol   713.00  ...    NaN  0.215149
999-21-3                diallyl maleate   693.00  ...  0.789  0.245054
999-52-0                3-chloroheptane   621.08  ...  0.418  0.248495
999-78-0         4,4-dimethyl-2-pentyne   552.81  ...  0.191  0.281303
999-97-3           hexamethyldisilazane   544.00  ...  0.510  0.260213
[7549 rows x 6 columns]

An extensive set of functions for calculating all sorts of chemical properties 
are available along with their respective coefficients for a wide range of chemicals:

>>> from chemicals import Antoine
>>> from chemicals.vapor_pressure import AntoinePoling
>>> antoine_coefficients = AntoinePoling.loc[CAS_water] # For calculating saturated vapor pressure 
>>> antoine_coefficients
Chemical     water 
A              10.1
B          1.69e+03
C               -43
Tmin            273
Tmax            473
Name: 7732-18-5, dtype: object
>>> _, A, B, C, *_ = antoine_coefficients
>>> T = 373.15 # Temperature [K]
>>> Antoine(T, A, B, C) # Vapor pressure [Pa]
101047.25357066597

Roadmap
-------

The authors' main development item is to provide the latest methods
for the calculation of chemical properties (both thermodynamic and transport) 
and extending the local databank.

# TODO: Do we want to mention phase equilibrium and flash algorithms?

# Should we move this to the Developer's guide?

Latest source code
------------------

The latest development version of chemicals's sources can be obtained at

    https://github.com/CalebBell/chemicals


Bug reports
-----------

To report bugs, please use the chemicals's Bug Tracker at:

    https://github.com/CalebBell/chemicals/issues


License information
-------------------

See ``LICENSE.txt`` for information on the terms & conditions for usage
of this software, and a DISCLAIMER OF ALL WARRANTIES.

Although not required by the chemicals license, if it is convenient for you,
please cite chemicals if used in your work. Please also consider contributing
any changes you make back, and benefit the community.


Citation
--------

# TODO: Do we want to keep this name and change thermo's description?

To cite chemicals in publications use:

    Caleb Bell and Contributors (2016-2020). chemicals: Chemical properties component of Chemical Engineering Design Library (ChEDL)
    https://github.com/CalebBell/chemicals.
