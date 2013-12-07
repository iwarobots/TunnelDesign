#!/usr/bin/env python


from __future__ import absolute_import, division

from math import radians

import nose
import nose.tools as nt

from properties.prandtl_meyer_function import (m2nu,
                                               nu2m)


@nt.raises(ValueError)
def test_mach_lesser_than_one():
    m = 0.1
    m2nu(m)


def test_normal_mach():
    m1 = 1.5
    nt.assert_almost_equal(m2nu(m1), radians(11.9052), places=4)

    m2 = 2.6
    nt.assert_almost_equal(m2nu(m2), radians(41.4147), places=4)


@nt.raises(ValueError)
def test_inverse_prandtl_meyer_function_when_nu_is_less_than_zero():
    nu = -1
    nu2m(nu)

def test_inverse_prandtl_meyer_function():
    nu1 = radians(11.9052)
    nt.assert_almost_equal(nu2m(nu1), 1.5, places=1)


if __name__ == '__main__':
    nose.main()