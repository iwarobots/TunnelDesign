#!/usr/bin/env python

"""Test Mach angle functions.
"""

from __future__ import absolute_import, division

import nose
import nose.tools as nt

from properties.prandtl_meyer_function import m2mu_in_deg


@nt.raises(ValueError)
def test_mach_lesser_than_one():
    m = 0.1
    m2mu_in_deg(m)


def test_normal_mach():
    m1 = 1.5
    nt.assert_almost_equal(m2mu_in_deg(m1),  41.8103148, places=3)

    m2 = 2.6
    nt.assert_almost_equal(m2mu_in_deg(m2),  22.6198649, places=3)


if __name__ == '__main__':
    nose.main()