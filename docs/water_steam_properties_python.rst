Computing Properties of Water and Steam in Python
=================================================

Water is a very special substance. It is abundant, cheap, hydrating, and great for many engineering applications. Whatever your modeling goal, there is a good change you will require properties of water at various conditions.

There is an international association, IAPWS, which publishes and coordinates some of the best research on the properties of water. There is a special equation of state just for water developed by them that very accurately computes the properties of water, called IAPWS-95. There is also a "shortcut" version called IAPWS-97 which is faster to solve but has reduced accuracy and various discontinuities.


There are quite a few implementations of IAPWS-95 and IAPWS-97 out there. Besides the many commercial implementations, the are the following excellent open source ones:

* `iapws <https://github.com/jjgomera/iapws>`_ by Juan José Gómez Romera, GPL3 licensed, containing IAPWS-95 and IAPWS-97 among other standards. Implemented in Python.
* `CoolProp <https://github.com/CoolProp/CoolProp>`_ by Ian Bell, MIT licensed and containing IAPWS-95 and IAPWS-97 along with their transport properties. Implemented in C++ with an excellent interface to Python among other languages.
* `freesteam <http://freesteam.sourceforge.net/>`_ by John Pye, GPL3 licensed, containing most of IAPWS-97 and the transport properties. Implemented in C.

There are many more, but these are the best developed libraries that can be used from Python. Water is so common and present in so many calculations that for many applications it is important to make it as fast as possible. IAPWS-95 is conventionally slow; properties are requested at a specified temperature `T` and pressure `P`, but the equation of state's input variables are temperature and density! A numerical solver must be used in this case to find the density which yields the specified pressure. This density-solution procedure is normally the slowest part, although computing some properties requires many derivatives that can be slow also.

A good conventional density solver will take ~10-30 μs on a modern computer. Only the CPU clockspeed really matters for this calculation time. It was discovered that with the use of `Common subexpression elimination <https://en.wikipedia.org/wiki/Common_subexpression_elimination>`_, the calculation could be speed up quite a lot. Additionally, if the IAPWS-95 density solution is initialized by the IAPWS-97 explicit calculation (applicable most of the time but not always), a few more iterations can be saved.

The net result of these optimizations is a greatly improved density solve time - normally 2.5-4 μs when running with PyPy or Numba. The con to this approach is that the code is nearly unreadable, and it would not be possible to update the coefficients without rewriting the implementation. As IAPWS-95 is a static model which will be the best one available for many years to come, this is an acceptable trade off.