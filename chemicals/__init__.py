# -*- coding: utf-8 -*-
'''Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
Copyright (C) 2016, 2017, 2018, 2019, 2020 Caleb Bell <Caleb.Andrew.Bell@gmail.com>

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

import os

from . import utils
from . import critical
from . import elements
from . import dipole
from . import dippr
from . import temperature
from . import miscdata
from . import environment
from . import refractivity
from . import solubility
from . import lennard_jones
from . import vapor_pressure
from . import virial
from . import phase_change
from . import triple
from . import acentric
from . import exceptions

from .utils import *
from .critical import *
from .elements import *
from .dipole import *
from .dippr import *
from .temperature import *
from .miscdata import *
from .environment import *
from .refractivity import *
from .solubility import *
from .lennard_jones import *
from .vapor_pressure import *
from .virial import *
from .phase_change import *
from .triple import *
from .acentric import *
from .exceptions import *

__all__ = ['utils', 'critical', 'elements', 'dipole', 'dippr', 'temperature', 
           'miscdata', 'environment', 'refractivity', 'solubility', 
           'lennard_jones', 'vapor_pressure', 'virial', 'phase_change',
           'triple', 'exceptions', 'acentric']

__all__.extend(exceptions.__all__)
__all__.extend(critical.__all__)
__all__.extend(utils.__all__)
__all__.extend(elements.__all__)
__all__.extend(dipole.__all__)
__all__.extend(dippr.__all__)
__all__.extend(temperature.__all__)
__all__.extend(miscdata.__all__)
__all__.extend(environment.__all__)
__all__.extend(refractivity.__all__)
__all__.extend(solubility.__all__)
__all__.extend(lennard_jones.__all__)
__all__.extend(vapor_pressure.__all__)
__all__.extend(virial.__all__)
__all__.extend(phase_change.__all__)
__all__.extend(triple.__all__)
__all__.extend(acentric.__all__)

chemicals_dir = os.path.dirname(__file__)

__version__ = '0.1.1'
