chemicals: Chemical properties component of Chemical Engineering Design Library (ChEDL)
=======================================================================================

# TODO: Logo here?

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
modeling chemical processes. In fact, the `Biorefinery Simulation and Techno-Economic Analysis Modules (BioSTEAM) <https://biosteam.readthedocs.io/en/latest/>`_ 
is reliant on the chemicals library for the simulation of unit operations.

If you need to know something about a chemical, give chemicals a try.

.. toctree::
   :maxdepth: 1
   :numbered:
   :caption: Tutorial
   
   Overview
   Installation
   tutorial/Getting_started
   tutorial/Numba support
   tutorial/Numpy array support
   tutorial/Pint Quantities support

.. toctree::
   :maxdepth: 2
   :caption: API
   
   chemicals.acentric
   chemicals.combustion
   chemicals.critical
   chemicals.dipole
   chemicals.dippr
   chemicals.elements
   chemicals.environment
   chemicals.exceptions
   chemicals.flash_basic
   chemicals.heat_capacity
   chemicals.identifiers
   chemicals.interface
   chemicals.lennard_jones
   chemicals.miscdata
   chemicals.numba
   chemicals.permittivity
   chemicals.phase_change
   chemicals.rachford_rice
   chemicals.reaction
   chemicals.refractivity
   chemicals.safety
   chemicals.solubility
   chemicals.temperature
   chemicals.thermal_conductivity
   chemicals.triple
   chemicals.utils
   chemicals.vapor_pressure
   chemicals.vectorized
   chemicals.virial
   chemicals.viscosity
   chemicals.volume


.. toctree::
   :maxdepth: 1
   :caption: Developer's guide
   
   developers

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

To cite chemicals in publications use::

    Caleb Bell and Contributors (2016-2020). chemicals: Chemical properties component of Chemical Engineering Design Library (ChEDL)
    https://github.com/CalebBell/chemicals.
    
Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

