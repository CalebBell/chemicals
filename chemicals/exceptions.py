# -*- coding: utf-8 -*-
'''Chemical Engineering Design Library (ChEDL). Utilities for process modeling.
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
SOFTWARE.'''

__all__ = ['TrivialSolutionError',
           'PhaseCountReducedError', 
           'PhaseExistenceImpossible', 
           'OverspeficiedError']

class OverspeficiedError(Exception):
    '''Generic error to raise when multiple values are given, when only
    one option should be given.
    '''


class TrivialSolutionError(Exception):
    '''Error raised SS converges to trivial solution
    '''
    def __init__(self, message, comp_difference, iterations, err):
        super().__init__(message)
        self.comp_difference = comp_difference
        self.iterations = iterations
        self.err = err


class PhaseCountReducedError(Exception):
    '''Error raised SS inner flash loop says all Ks are under 1 or above 1
    '''
    def __init__(self, message, zs=None, Ks=None):
        super().__init__(message)
        self.zs = zs
        self.Ks = Ks


class PhaseExistenceImpossible(Exception):
    '''Error raised SS inner flash loop says all Ks are under 1 or above 1
    '''
    def __init__(self, message, zs=None, T=None, P=None):
        super().__init__(message)
        self.zs = zs
        self.T = T
        self.P = P
