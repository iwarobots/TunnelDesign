#!/usr/bin/env python


from __future__ import absolute_import, division

from math import asin, atan, degrees, radians, sqrt

from scipy.optimize import brentq

from properties.constants import GAMMA


MIN_MACH = 1E-5
MAX_MACH = 1E1


def m2nu_in_rad(m):
    if m < 1:
        raise ValueError('Mach number should be greater than or equal to 1')
    a = (GAMMA+1) / (GAMMA-1)
    b = m**2 - 1
    c = a**-1 * b
    return sqrt(a) * atan(sqrt(c)) - atan(sqrt(b))


def m2nu_in_deg(m):
    return degrees(m2nu_in_rad(m))


def nu_in_rad2m(nu):
    if nu < 0:
        raise ValueError('nu should be greater than or equal to 0')
    return brentq(lambda x: m2nu_in_rad(x)-nu, 1, MAX_MACH)


def nu_in_deg2m(nu):
    nu = radians(nu)
    return nu_in_rad2m(nu)


def m2mu_in_rad(m):
    return asin(1/m)


def m2mu_in_deg(m):
    return degrees(m2mu_in_rad(m))