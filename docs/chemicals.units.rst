Support for pint Quantities (chemicals.units)
=============================================

Basic module which wraps all chemicals functions and classes to be compatible with the
`pint <https://github.com/hgrecco/pint>`_ unit handling library.
All other object - dicts, lists, etc - are not wrapped. 

>>> import chemicals
>>> chemicals.units.Antoine # doctest: +ELLIPSIS
<function Antoine at 0x...>

The chemicals.units module also supports star imports; the same objects exported when importing from the main library
will be imported from chemicals.units.

>>> from chemicals.units import *

>>> CAS = CAS_from_any('methanol')
>>> Tc(CAS), Pc(CAS), Vc(CAS), Zc(CAS)
(<Quantity(512.5, 'kelvin')>, <Quantity(8084000.0, 'pascal')>, <Quantity(0.000117, 'meter ** 3 / mole')>, <Quantity(0.223, 'dimensionless')>)
>>> (Tt(CAS), Tm(CAS), Tb(CAS), Pt(CAS))
(<Quantity(175.59, 'kelvin')>, <Quantity(175.15, 'kelvin')>, <Quantity(337.65, 'kelvin')>, None)

>>> iapws95_rho(T=55*u.degF, P=500*u.psi)
<Quantity(1000.97992, 'kilogram / meter ** 3')>
>>> sigma_IAPWS(200*u.degR)
<Quantity(0.0897667127, 'newton / meter')>

>>> molecular_weight({'H': 12, 'C': 20, 'O': 5})
<Quantity(332.30628, 'gram / mole')>

Functions that do not return numbers are not converted into pint quantities, for example:

>>> atoms_to_Hill({'H': 5, 'C': 2, 'Br': 1})
'C2H5Br'

Functions that return dimensionless numbers are pint quantities.

>>> logP('67-56-1')
<Quantity(-0.74, 'dimensionless')>
>>> Stiel_polar_factor(Psat=169745*u.Pa, Pc=22048321.0*u.Pa, omega=0.344)
<Quantity(0.0232214674, 'dimensionless')>

It is also possible to use `chemicals.units` without the star import:

>>> import chemicals.units

When a function is used with inputs that should have units but they aren't provided by the user, an error is raised.

>>> ideal_gas(298.15, 101325.) # doctest: +IGNORE_EXCEPTION_DETAIL
Traceback (most recent call last):
TypeError: 298.15 has no quantity


For further information on this interface, please see the documentation of `fluids.units <https://fluids.readthedocs.io/fluids.units.html>`_ which is built in the same way.