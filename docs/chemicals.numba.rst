Support for Numba (chemicals.numba)
===================================

Basic module which wraps most of chemicals functions and classes to be compatible with the
`Numba <https://github.com/numba/numba>`_ dynamic Python compiler.
Numba is only supported on Python 3, and may require the latest version of Numba.
Numba is rapidly evolving, and hopefully in the future it will support more of
the functionality of fluids.

Using the numba-accelerated version of `chemicals` is easy; simply call functions
and classes from the fluids.numba namespace. The chemicals.numba module must be
imported separately; it is not loaded automatically as part of fluids.

>>> import chemicals
>>> import chemicals.numba
>>> chemicals.numba.Antoine(180, A=8.95894, B=510.595, C=-15.95)
702271.0518579542

There is a delay while the code is compiled when using Numba;
the speed is not quite free.


Some parts of `chemicals` are not supported in the Numba interface.
Any data lookup function like :py:func:`~.Tc`  and :py:func:`~.omega_methods` are not numba compatible.
These functions could not be speed up anyway though as they are memory-bound
and numba helps with things that are compute-bound.

Today, the list of things known not to work is as follows:

- Everything in :py:mod:`chemicals.identifiers`
- Everything in :py:mod:`chemicals.combustion` (uses dictionaries)
- Everything in :py:mod:`chemicals.elements` (uses dictionaries, periodic_table class, string parsing)
- Everything in :py:mod:`chemicals.critical` except the mixture critical point routines
- In :py:mod:`chemicals.flash_basic`, :py:func:`~.flash_Tb_Tc_Pc` and :py:func:`~.flash_ideal` have not been ported but can be made compatible.
- In :py:mod:`chemicals.dippr`, only :py:func:`~.EQ102` is unsupported as it uses a complex hyp2f1 call.

- The assorted functions  :py:func:`~.TRCCp`  :py:func:`~.Lastovka_Shaw_T_for_Sm`  :py:func:`~.Lastovka_Shaw_T_for_Hm`  :py:func:`~.collision_integral_Neufeld_Janzen_Aziz`  :py:func:`~.collision_integral_Kim_Monroe`  :py:func:`~.Gibbs_formation`  :py:func:`~.stoichiometric_matrix`   :py:func:`~.balance_stoichiometry`   :py:func:`~.NFPA_combustible_classification`  :py:func:`~.fire_mixing`   :py:func:`~.Henry_pressure_mixture`   :py:func:`~.Henry_converter`   :py:func:`~.T_converter`   :py:func:`~.Z_from_virial_density_form`  :py:func:`~.Z_from_virial_pressure_form`  :py:func:`~.Wilke`  :py:func:`~.Wilke_prefactors` :py:func:`~.Brokaw` :py:func:`~.viscosity_converter` :py:func:`~.Wilke_prefactors` :py:func:`~.Wilke_prefactors`  may or may not be able to be ported but have not yet been.


All of the regular Numba-compiled functions are built with the `nogil` flag,
which means you can use Python's threading mechanism effectively to get
the speed of parallel processing.
