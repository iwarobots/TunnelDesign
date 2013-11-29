#!/usr/bin/env python


from math import atan, pi, sqrt

from properties.constants import GAMMA


def nu_in_rad(m):
    if m < 1:
        raise ValueError('Mach number should be greater than or equal to 1')
    a = (GAMMA+1) / (GAMMA-1)
    b = m**2 - 1
    c = a**-1 * b
    return sqrt(a) * atan(sqrt(c)) - atan(sqrt(b))


def nu_in_deg(m):
    return nu_in_rad(m) * 180 / pi
