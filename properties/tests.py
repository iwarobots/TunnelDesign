#!/usr/bin/env python


import nose
import nose.tools as nt

from properties.prandtl_meyer_function import nu_in_deg


@nt.raises(ValueError)
def test_mach_lesser_than_one():
    m = 0.1
    nu_in_deg(m)


def test_normal_mach():
    m1 = 1.5
    nt.assert_almost_equal(nu_in_deg(m1), 11.9052, places=4)

    m2 = 2.6
    nt.assert_almost_equal(nu_in_deg(m2), 41.4147, places=4)


if __name__ == '__main__':
    nose.main()