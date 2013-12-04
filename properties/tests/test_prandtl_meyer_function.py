#!/usr/bin/env python


from __future__ import absolute_import, division

import nose
import nose.tools as nt

from properties.prandtl_meyer_function import (m2nu_in_deg,
                                               nu_in_deg2m)


@nt.raises(ValueError)
def test_mach_lesser_than_one():
    m = 0.1
    m2nu_in_deg(m)


def test_normal_mach():
    m1 = 1.5
    nt.assert_almost_equal(m2nu_in_deg(m1), 11.9052, places=4)

    m2 = 2.6
    nt.assert_almost_equal(m2nu_in_deg(m2), 41.4147, places=4)


@nt.raises(ValueError)
def test_inverse_prandtl_meyer_function_when_nu_is_less_than_zero():
    nu = -1
    nu_in_deg2m(nu)

def test_inverse_prandtl_meyer_function():
    nu1 = 11.9052
    nt.assert_almost_equal(nu_in_deg2m(nu1), 1.5, places=1)


if __name__ == '__main__':
    nose.main()