from fluids.constants import *

#IS_PYPY = True
from fluids.numerics import IS_PYPY, normalize

if not IS_PYPY:
    import chemicals.numba
import numpy as np


def also_numba(f):
    if not IS_PYPY:
        f.duplicate_with_numba = True
    return f



class BaseTimeSuite:
    def setup(self):
        if not IS_PYPY:
            for k in dir(self.__class__):
                if 'time' in k and 'numba' in k:
                    c = getattr(self, k)
                    c()

from chemicals.viscosity import *

if not IS_PYPY:
    mu_IAPWS_numba = chemicals.numba.mu_IAPWS
    mu_air_lemmon_numba = chemicals.numba.mu_air_lemmon
    PPDS9_numba = chemicals.numba.PPDS9
    dPPDS9_dT_numba = chemicals.numba.dPPDS9_dT
    Letsou_Stiel_numba = chemicals.numba.Letsou_Stiel
    Przedziecki_Sridhar_numba = chemicals.numba.Przedziecki_Sridhar
    Lucas_numba = chemicals.numba.Lucas
    Yoon_Thodos_numba = chemicals.numba.Yoon_Thodos
    Stiel_Thodos_numba = chemicals.numba.Stiel_Thodos
    Lucas_gas_numba = chemicals.numba.Lucas_gas
    viscosity_gas_Gharagheizi_numba = chemicals.numba.viscosity_gas_Gharagheizi
    Twu_1985_numba = chemicals.numba.Twu_1985
    viscosity_index_numba = chemicals.numba.viscosity_index

class TimeViscositySuite(BaseTimeSuite):
    def time_mu_IAPWS(self):
        mu_IAPWS(298.15, 998.)
    def time_mu_IAPWS_numba(self):
        mu_IAPWS_numba(298.15, 998.)
    def time_mu_IAPWS_full(self):
        mu_IAPWS(T=647.35, rho=222, drho_dP=175.456980972231e-6, drho_dP_Tr=3.119177410324e-6)
    def time_mu_IAPWS_full_numba(self):
        mu_IAPWS_numba(T=647.35, rho=222, drho_dP=175.456980972231e-6, drho_dP_Tr=3.119177410324e-6)

    def time_mu_air_lemmon(self):
        mu_air_lemmon(298.15, 40.)
    def time_mu_air_lemmon_numba(self):
        mu_air_lemmon_numba(298.15, 40.)

    def time_PPDS9(self):
        PPDS9(400.0, 1.74793, 1.33728, 482.347, 41.78, 9.963e-05)
    def time_PPDS9_numba(self):
        PPDS9_numba(400.0, 1.74793, 1.33728, 482.347, 41.78, 9.963e-05)
    def time_dPPDS9_dT(self):
        dPPDS9_dT(400.0, 1.74793, 1.33728, 482.347, 41.78, 9.963e-05)
    def time_dPPDS9_dT_numba(self):
        dPPDS9_dT_numba(400.0, 1.74793, 1.33728, 482.347, 41.78, 9.963e-05)

    def time_Letsou_Stiel(self):
        Letsou_Stiel(400., 46.07, 516.25, 6.383E6, 0.6371)
    def time_Letsou_Stiel_numba(self):
        Letsou_Stiel_numba(400., 46.07, 516.25, 6.383E6, 0.6371)

    def time_Przedziecki_Sridhar(self):
        Przedziecki_Sridhar(383., 178., 591.8, 41E5, 316E-6, 95E-6, .263, 92.14)
    def time_Przedziecki_Sridhar_numba(self):
        Przedziecki_Sridhar_numba(383., 178., 591.8, 41E5, 316E-6, 95E-6, .263, 92.14)

    def time_Lucas(self):
        Lucas(300., 500E5, 572.2, 34.7E5, 0.236, 0, 0.00068)
    def time_Lucas_numba(self):
        Lucas_numba(300., 500E5, 572.2, 34.7E5, 0.236, 0, 0.00068)

    def time_Yoon_Thodos(self):
        Yoon_Thodos(300., 556.35, 4.5596E6, 153.8)
    def time_Yoon_Thodos_numba(self):
        Yoon_Thodos_numba(300., 556.35, 4.5596E6, 153.8)

    def time_Stiel_Thodos(self):
        Stiel_Thodos(300., 556.35, 4.5596E6, 153.8)
    def time_Stiel_Thodos_numba(self):
        Stiel_Thodos_numba(300., 556.35, 4.5596E6, 153.8)

    def time_Lucas_gas(self):
        Lucas_gas(T=550., Tc=512.6, Pc=80.9E5, Zc=0.224, MW=32.042, dipole=1.7)
    def time_Lucas_gas_numba(self):
        Lucas_gas_numba(T=550., Tc=512.6, Pc=80.9E5, Zc=0.224, MW=32.042, dipole=1.7)

    def time_viscosity_gas_Gharagheizi(self):
        viscosity_gas_Gharagheizi(120., 190.564, 45.99E5, 16.04246)
    def time_viscosity_gas_Gharagheizi_numba(self):
        viscosity_gas_Gharagheizi_numba(120., 190.564, 45.99E5, 16.04246)

    def time_Twu_1985(self):
        Twu_1985(T=338.7055, Tb=672.3166, rho=895.5189)
    def time_Twu_1985_numba(self):
        Twu_1985_numba(T=338.7055, Tb=672.3166, rho=895.5189)

    def time_viscosity_index(self):
        viscosity_index(73.3E-6, 8.86E-6, rounding=True)
    def time_viscosity_index_numba(self):
        viscosity_index_numba(73.3E-6, 8.86E-6, rounding=True)

    def time_viscosity_converter_1(self):
        viscosity_converter(8.79, 'engler', 'parlin cup #7')
    def time_viscosity_converter_2(self):
        viscosity_converter(700, 'Saybolt Universal Seconds', 'kinematic viscosity')

from chemicals.permittivity import *

if not IS_PYPY:
    permittivity_IAPWS_numba = chemicals.numba.permittivity_IAPWS

class TimePermittivitySuite(BaseTimeSuite):
    def time_permittivity_IAPWS(self):
        permittivity_IAPWS(650., 40.31090)
    def time_permittivity_IAPWS_numba(self):
        permittivity_IAPWS_numba(650., 40.31090)


from chemicals.thermal_conductivity import *

if not IS_PYPY:
    k_IAPWS_numba = chemicals.numba.k_IAPWS
    k_air_lemmon_numba = chemicals.numba.k_air_lemmon
    Sheffy_Johnson_numba = chemicals.numba.Sheffy_Johnson
    Sato_Riedel_numba = chemicals.numba.Sato_Riedel
    Gharagheizi_liquid_numba = chemicals.numba.Gharagheizi_liquid
    Nicola_original_numba = chemicals.numba.Nicola_original
    Lakshmi_Prasad_numba = chemicals.Lakshmi_Prasad
    Nicola_numba = chemicals.numba.Nicola
    Bahadori_liquid_numba = chemicals.numba.Bahadori_liquid
    kl_Mersmann_Kind_numba = chemicals.numba.kl_Mersmann_Kind
    DIPPR9G_numba = chemicals.numba.DIPPR9G
    Missenard_numba = chemicals.numba.Missenard
    Eucken_numba = chemicals.numba.Eucken
    Eucken_modified_numba = chemicals.numba.Eucken_modified
    DIPPR9B_numba = chemicals.numba.DIPPR9B
    Chung_numba = chemicals.numba.Chung



class TimeThermalConductivitySuite(BaseTimeSuite):
    def time_k_IAPWS(self):
        k_IAPWS(647.35, 750.)
    def time_k_IAPWS_numba(self):
        k_IAPWS_numba(647.35, 750.)
    def time_k_IAPWS_full(self):
        k_IAPWS(T=620., rho=613.227777440324, Cp=7634.337046792, Cv=3037.934412104, mu=70.905106751524E-6, drho_dP=5.209378197916E-6)
    def time_k_IAPWS_full_numba(self):
        k_IAPWS_numba(T=620., rho=613.227777440324, Cp=7634.337046792, Cv=3037.934412104, mu=70.905106751524E-6, drho_dP=5.209378197916E-6)

    def time_k_air_lemmon(self):
        k_air_lemmon(300.0, 40.0)
    def time_k_air_lemmon_numba(self):
        k_air_lemmon_numba(300.0, 40.0)
    def time_k_air_lemmon_full(self):
        k_air_lemmon(132.64, 10400, 2137.078854678728, 35.24316159996235, 0.07417878614315769, 0.00035919027241528256, 1.7762253265868595e-05)
    def time_k_air_lemmon_full_numba(self):
        k_air_lemmon_numba(132.64, 10400, 2137.078854678728, 35.24316159996235, 0.07417878614315769, 0.00035919027241528256, 1.7762253265868595e-05)

    def time_Sheffy_Johnson(self):
        Sheffy_Johnson(300, 47, 280)
    def time_Sheffy_Johnson_numba(self):
        Sheffy_Johnson_numba(300, 47, 280)

    def time_Sato_Riedel(self):
        Sato_Riedel(300, 47, 390, 520)
    def time_Sato_Riedel_numba(self):
        Sato_Riedel_numba(300, 47, 390, 520)

    def time_Lakshmi_Prasad(self):
        Lakshmi_Prasad(273.15, 100)
    def time_Lakshmi_Prasad_numba(self):
        Lakshmi_Prasad_numba(273.15, 100)

    def time_Gharagheizi_liquid(self):
        Gharagheizi_liquid(300, 40, 350, 1E6, 0.27)
    def time_Gharagheizi_liquid_numba(self):
        Gharagheizi_liquid_numba(300, 40, 350, 1E6, 0.27)

    def time_Nicola_original(self):
        Nicola_original(300, 142.3, 611.7, 0.49, 201853)
    def time_Nicola_original_numba(self):
        Nicola_original_numba(300, 142.3, 611.7, 0.49, 201853)

    def time_Nicola(self):
        Nicola(300, 142.3, 611.7, 0.49, 201853)
    def time_Nicola_numba(self):
        Nicola_numba(300, 142.3, 611.7, 0.49, 201853)

    def time_Bahadori_liquid(self):
        Bahadori_liquid(273.15, 170)
    def time_Bahadori_liquid_numba(self):
        Bahadori_liquid_numba(273.15, 170)

    def time_kl_Mersmann_Kind(self):
        kl_Mersmann_Kind(400, 170.33484, 658.0, 0.000754, 38)
    def time_kl_Mersmann_Kind_numba(self):
        kl_Mersmann_Kind_numba(400, 170.33484, 658.0, 0.000754, 38)

    def time_DIPPR9G(self):
        DIPPR9G(515.05, 3.92E7, 579.15, 3.212E6, 7.085E-2)
    def time_DIPPR9G_numba(self):
        DIPPR9G_numba(515.05, 3.92E7, 579.15, 3.212E6, 7.085E-2)

    def time_Missenard(self):
        Missenard(304., 6330E5, 591.8, 41E5, 0.129)
    def time_Missenard_numba(self):
        Missenard_numba(304., 6330E5, 591.8, 41E5, 0.129)

    def time_Eucken(self):
        Eucken(MW=72.151, Cvm=135.9, mu=8.77E-6)
    def time_Eucken_numba(self):
        Eucken_numba(MW=72.151, Cvm=135.9, mu=8.77E-6)

    def time_Eucken_modified(self):
        Eucken_modified(MW=72.151, Cvm=135.9, mu=8.77E-6)
    def time_Eucken_modified_numba(self):
        Eucken_modified_numba(MW=72.151, Cvm=135.9, mu=8.77E-6)

    def time_DIPPR9B(self):
        DIPPR9B(200., 28.01, 20.826, 1.277E-5, 132.92, chemtype='linear')
    def time_DIPPR9B_numba(self):
        DIPPR9B_numba(200., 28.01, 20.826, 1.277E-5, 132.92, chemtype='linear')


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
    Weinaug_Katz_numba = chemicals.numba.Weinaug_Katz
    Winterfeld_Scriven_Davis_numba = chemicals.numba.Winterfeld_Scriven_Davis
    Diguilio_Teja_numba = chemicals.numba.Diguilio_Teja


Diguilio_Teja_ns = Winterfeld_Scriven_Davis_ns = Weinaug_Katz_ns = (2, 5, 10, 20, 50, 100)
class TimeInterfaceSuite(BaseTimeSuite):
    def __init__(self):
        super().__init__()

        for N in Diguilio_Teja_ns:
            Diguilio_Teja_kwargs = dict(T=298.15, xs=normalize([0.1606, 0.8394]*N), sigmas_Tb=[0.01424, 0.02530]*N, Tbs=[309.21, 312.95]*N, Tcs=[469.7, 508.0]*N)
            N *= 2
            setattr(self, 'DT%d' %N, Diguilio_Teja_kwargs)
            Diguilio_Teja_kwargs = Diguilio_Teja_kwargs.copy()
            for s in ('sigmas_Tb', 'xs', 'Tbs', 'Tcs'):
                Diguilio_Teja_kwargs[s] = np.array(Diguilio_Teja_kwargs[s])
            setattr(self, 'DTnp%d' %N, Diguilio_Teja_kwargs)


        for N in Weinaug_Katz_ns:
            Weinaug_Katz_kwargs = dict(parachors=[5.1e-5, 7.2e-5]*N, Vml=0.000125, Vmg=0.02011, xs=normalize([.4, .6]*N), ys=normalize([.6, .4]*N))
            N *= 2
            setattr(self, 'WK%d' %N, Weinaug_Katz_kwargs)
            Weinaug_Katz_kwargs = Weinaug_Katz_kwargs.copy()
            for s in ('parachors', 'xs', 'ys'):
                Weinaug_Katz_kwargs[s] = np.array(Weinaug_Katz_kwargs[s])
            setattr(self, 'WKnp%d' %N, Weinaug_Katz_kwargs)

        for N in Winterfeld_Scriven_Davis_ns:
            Winterfeld_Scriven_Davis_kwargs = dict(xs=normalize([0.1606, 0.8394]*N), sigmas=[0.01547, 0.02877]*N, rhoms=[8610., 15530.]*N)
            N *= 2
            setattr(self, 'WS%d' %N, Winterfeld_Scriven_Davis_kwargs)
            Winterfeld_Scriven_Davis_kwargs = Winterfeld_Scriven_Davis_kwargs.copy()
            for s in ('sigmas', 'xs', 'rhoms'):
                Winterfeld_Scriven_Davis_kwargs[s] = np.array(Winterfeld_Scriven_Davis_kwargs[s])
            setattr(self, 'WSnp%d' %N, Winterfeld_Scriven_Davis_kwargs)


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

for n in Weinaug_Katz_ns:
    n *= 2
    string = 'WK%d' %(n,)
    def f(self, N=n, string=string):
        kwargs = getattr(self, string)
        Weinaug_Katz(**kwargs)
    setattr(TimeInterfaceSuite, 'time_Weinaug_Katz_%d' %(n,), f)

    string = 'WKnp%d' %(n,)
    def fnp(self, N=n, string=string):
        kwargs = getattr(self, string)
        Weinaug_Katz_numba(**kwargs)
    setattr(TimeInterfaceSuite, 'time_Weinaug_Katz_%d_numba' %(n,), fnp)

for n in Winterfeld_Scriven_Davis_ns:
    n *= 2
    string = 'WS%d' %(n,)
    def f(self, N=n, string=string):
        kwargs = getattr(self, string)
        Winterfeld_Scriven_Davis(**kwargs)
    setattr(TimeInterfaceSuite, 'time_Winterfeld_Scriven_Davis_%d' %(n,), f)

    string = 'WSnp%d' %(n,)
    def fnp(self, N=n, string=string):
        kwargs = getattr(self, string)
        Winterfeld_Scriven_Davis_numba(**kwargs)
    setattr(TimeInterfaceSuite, 'time_Winterfeld_Scriven_Davis_%d_numba' %(n,), fnp)

for n in Diguilio_Teja_ns:
    n *= 2
    string = 'DT%d' %(n,)
    def f(self, N=n, string=string):
        kwargs = getattr(self, string)
        Diguilio_Teja(**kwargs)
    setattr(TimeInterfaceSuite, 'time_Diguilio_Teja_%d' %(n,), f)

    string = 'DTnp%d' %(n,)
    def fnp(self, N=n, string=string):
        kwargs = getattr(self, string)
        Diguilio_Teja_numba(**kwargs)
    setattr(TimeInterfaceSuite, 'time_Diguilio_Teja_%d_numba' %(n,), fnp)


from chemicals.dippr import *

if IS_PYPY:
    EQ100_numba = chemicals.numba.EQ100
    EQ101_numba = chemicals.numba.EQ101
    EQ102_numba = chemicals.numba.EQ102
    EQ104_numba = chemicals.numba.EQ104
    EQ105_numba = chemicals.numba.EQ105
    EQ106_numba = chemicals.numba.EQ106
    EQ107_numba = chemicals.numba.EQ107
    EQ114_numba = chemicals.numba.EQ114
    EQ115_numba = chemicals.numba.EQ115
    EQ116_numba = chemicals.numba.EQ116
    EQ127_numba = chemicals.numba.EQ127

suites = [TimeInterfaceSuite, TimePermittivitySuite, TimeViscositySuite,
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

