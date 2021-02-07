from fluids.numerics import IS_PYPY
from fluids.constants import *
#IS_PYPY = True

if not IS_PYPY:
    import fluids.numba
    import chemicals.numba
    import numba
from datetime import datetime
import numpy as np

def also_numba(f):
    if not IS_PYPY:
        f.duplicate_with_numba = True
    return f



class BaseTimeSuite(object):
    def setup(self):
        if not IS_PYPY:
            for k in dir(self.__class__):
                if 'time' in k and 'numba' in k:
                    c = getattr(self, k)
                    c()


from chemicals.interface import * 

if not IS_PYPY:
    sigma_IAPWS_numba = chemicals.numba.sigma_IAPWS
    REFPROP_sigma_numba = chemicals.numba.REFPROP_sigma
    Somayajulu_numba = chemicals.numba.Somayajulu
    Brock_Bird_numba = chemicals.numba.Brock_Bird
    Pitzer_sigma_numba = chemicals.numba.Pitzer_sigma
    Sastri_Rao_numba = chemicals.numba.Sastri_Rao
    Zuo_Stenby_numba = chemicals.numba.Zuo_Stenby
    Hakim_Steinberg_Stiel_numba = chemicals.numba.Hakim_Steinberg_Stiel
    Miqueu_numba = chemicals.numba.Miqueu
    Aleem_numba = chemicals.numba.Aleem
    Mersmann_Kind_sigma_numba = chemicals.numba.Mersmann_Kind_sigma
    API10A32_numba = chemicals.numba.API10A32
    Meybodi_Daryasafar_Karimi_numba = chemicals.numba.Meybodi_Daryasafar_Karimi

class TimeInterfaceSuite(BaseTimeSuite):
    def time_sigma_IAPWS(self):
        sigma_IAPWS(450.)
    def time_sigma_IAPWS_numba(self):
        sigma_IAPWS_numba(450.)
        
        
    def time_REFPROP_sigma(self):
        REFPROP_sigma(298.15, 647.096, -0.1306, 2.471, 0.2151, 1.233)
    def time_REFPROP_sigma_numba(self):
        REFPROP_sigma_numba(298.15, 647.096, -0.1306, 2.471, 0.2151, 1.233)
        
        
    def time_Somayajulu(self):
        Somayajulu(300, 647.126, 232.713514, -140.18645, -4.890098)
    def time_Somayajulu_numba(self):
        Somayajulu_numba(300, 647.126, 232.713514, -140.18645, -4.890098)
        
    def time_Brock_Bird(self):
        Brock_Bird(293.15, 404.75, 633.0, 4530000.0)
    def time_Brock_Bird_numba(self):
        Brock_Bird_numba(293.15, 404.75, 633.0, 4530000.0)
        
    def time_Pitzer_sigma(self):
        Pitzer_sigma(293., 633.0, 4530000.0, 0.249)
    def time_Pitzer_sigma_numba(self):
        Pitzer_sigma_numba(293., 633.0, 4530000.0, 0.249)
        
    def time_Sastri_Rao(self):
        Sastri_Rao(293.15, 404.75, 633.0, 4530000.0)
    def time_Sastri_Rao_numba(self):
        Sastri_Rao_numba(293.15, 404.75, 633.0, 4530000.0)
        
    def time_Zuo_Stenby(self):
        Zuo_Stenby(293., 633.0, 4530000.0, 0.249)
    def time_Zuo_Stenby_numba(self):
        Zuo_Stenby_numba(293., 633.0, 4530000.0, 0.249)
        
    def time_Hakim_Steinberg_Stiel(self):
        Hakim_Steinberg_Stiel(298.15, 563.0, 4414000.0, 0.59, StielPolar=-0.07872)
    def time_Hakim_Steinberg_Stiel_numba(self):
        Hakim_Steinberg_Stiel_numba(298.15, 563.0, 4414000.0, 0.59, StielPolar=-0.07872)
        
    def time_Miqueu(self):
        Miqueu(300., 340.1, 0.000199, 0.1687)
    def time_Miqueu_numba(self):
        Miqueu_numba(300., 340.1, 0.000199, 0.1687)
        
    def time_Aleem(self):
        Aleem(T=90, MW=16.04246, Tb=111.6, rhol=458.7, Hvap_Tb=510870., Cpl=2465.)
    def time_Aleem_numba(self):
        Aleem_numba(T=90, MW=16.04246, Tb=111.6, rhol=458.7, Hvap_Tb=510870., Cpl=2465.)
        
    def time_Mersmann_Kind_sigma(self):
        Mersmann_Kind_sigma(298.15, 164.15, 328.25, 497.1, 3430000.0)
    def time_Mersmann_Kind_sigma_numba(self):
        Mersmann_Kind_sigma_numba(298.15, 164.15, 328.25, 497.1, 3430000.0)
        
    def time_API10A32(self):
        API10A32(288.7, 741.4, 12.4)
    def time_API10A32_numba(self):
        API10A32_numba(288.7, 741.4, 12.4)
        
    def time_Meybodi_Daryasafar_Karimi(self):
        Meybodi_Daryasafar_Karimi(980, 760, 580, 914)
    def time_Meybodi_Daryasafar_Karimi_numba(self):
        Meybodi_Daryasafar_Karimi_numba(980, 760, 580, 914)

suites = [TimeInterfaceSuite,
          ]


                
for suite in suites:
    continue
    # asv requires inspect to work :(
    # Do I want to write a file that writes this benchmark file?
    glbs, lcls = {}, {}
    for k in dir(suite):
        if 'time' in k:
            f = getattr(suite, k)
            if hasattr(f, 'duplicate_with_numba'):
                source = inspect.getsource(f)
                source = '\n'.join([s[4:] for s in source.split('\n')[1:]])
                orig_function = k.replace('time_', '')
                numba_function = orig_function + '_numba'
                new_function_name = k + '_numba'
                new_source = source.replace(orig_function, numba_function)
                exec(new_source, glbs, lcls)
                setattr(suite, new_function_name, lcls[new_function_name])

if IS_PYPY:
    for s in suites:
        for k in dir(s):
            if 'time' in k and 'numba' in k:
                delattr(s, k)
                
