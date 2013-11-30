#!/usr/bin/env python

"""Test Mach angle functions.

Test data is obtained from http://www.grc.nasa.gov/WWW/k-12/airplane/machang.html.
"""


import nose
import nose.tools as nt

from properties.prandtl_meyer_function import mu_in_deg


@nt.raises(ValueError)
def test_mach_lesser_than_one():
    m = 0.1
    mu_in_deg(m)


def test_normal_mach():
    m1 = 1.5
    nt.assert_almost_equal(mu_in_deg(m1), 41.762, places=3)

    m2 = 2.6
    nt.assert_almost_equal(mu_in_deg(m2), 22.594, places=3)


if __name__ == '__main__':
    nose.main()