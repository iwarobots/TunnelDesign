#!/usr/bin/env python

# TODO: Needs tests
# TODO: Needs performance improvements


from __future__ import absolute_import, division

from math import atan, tan

import numpy as np

from scipy.optimize import brentq
from scipy.linalg import solve

import matplotlib.pyplot as plt

import properties.prandtl_meyer_function as pmf


def angle_between_two_lines(k1, k2):
    return atan(abs((k2-k1)/(1+k1*k2)))


def intersection_of_two_rays(x1, y1, k1, x2, y2, k2):
    a = np.array([[k1, -1],
                  [k2, -1]])
    b = np.array([k1*x1-y1, k2*x2-y2])
    return solve(a, b)


class Solver(object):

    def __init__(self,
                 final_mach,
                 n):

        # Checks invalid input params
        # Checks Mach number
        if final_mach <= 1:
            raise ValueError('Mach number at exit should be greater than 1')
        
        # Checks expansion steps
        if n <= 1:
            raise ValueError('Expansion steps should be greater than 1')

        self._final_mach = final_mach
        self._n = n

        self._size = n + 1
        self._region_mat = None
        self._wall_regions = None
        self._points = None
        self._theta_plus_nu = None
        self._theta_minus_nu = None
        self._theta = None
        self._nu = None
        self._m = None
        self._mu = None

    @property
    def final_mach(self):
        return self._final_mach

    @property
    def n(self):
        return self._n

    @property
    def size(self):
        return self._size

    @property
    def region_mat(self):
        if self._region_mat is None:
            self._region_mat = self.compute_region_mat()
        return self._region_mat

    @property
    def wall_regions(self):
        if self._wall_regions is None:
            self._wall_regions = self.compute_wall_regions()
        return self._wall_regions

    @property
    def theta_max(self):
        return pmf.m2nu_in_rad(self.final_mach) / 2

    @property
    def theta_per_step(self):
        return self.theta_max / self.n

    def compute_region_mat(self):
        self._region_mat = np.zeros((self.size, self.size), dtype=int)
        self._region_mat[:] = np.nan

        cnt = 0
        for i in xrange(self.size):
            for j in xrange(i, self.size):
                self._region_mat[i, j] = cnt
                cnt += 1
        return self._region_mat

    def compute_wall_regions(self):
        return self.region_mat[:, self.n]

    def compute_theta_minus_nu(self):
        self._theta_minus_nu = np.zeros((self.size, self.size))
        self._theta_minus_nu[:] = np.nan

        for i in xrange(self.size):
            for j in xrange(i, self.size):
                self._theta_minus_nu[i, j] = -2 * i * self.theta_per_step
        return self._theta_minus_nu

    def compute_theta_plus_nu(self):
        self._theta_plus_nu = np.zeros((self.size, self.size))
        self._theta_plus_nu[:] = np.nan

        for i in xrange(self.size):
            for j in xrange(i+1):
                self._theta_plus_nu[i, j] = 2 * i * self.theta_per_step
        self._theta_plus_nu = self._theta_plus_nu.T
        return self._theta_plus_nu

    def _tpn_tmn_is_computed(self):
        if self._theta_plus_nu is None:
            self.compute_theta_plus_nu()
        if self._theta_minus_nu is None:
            self.compute_theta_minus_nu()

    def compute_theta(self):
        self._tpn_tmn_is_computed()
        self._theta = (self._theta_plus_nu+self._theta_minus_nu) / 2
        return self._theta

    def compute_nu(self):
        self._tpn_tmn_is_computed()
        self._nu = (self._theta_plus_nu-self._theta_minus_nu) / 2
        return self._nu

    def compute_m(self):
        if self._nu is None:
            self.compute_nu()

        self._m = np.zeros((self.size, self.size))
        self._m[:] = np.nan

        for i in xrange(self.size):
            for j in xrange(i, self.size):
                self._m[i, j] = pmf.nu_in_rad2m(self._nu[i, j])
        return self._m

    def compute_mu(self):
        if self._m is None:
            self.compute_m()

        self._mu = np.zeros((self.size, self.size))
        self._mu[:] = np.nan

        for i in xrange(self.size):
            for j in xrange(i, self.size):
                self._mu[i, j] = pmf.m2mu_in_rad(self._m[i, j])
        return self._mu

    def _compute_axis_point(self, n):
        if not 0 <= n <= self._n - 1:
            raise ValueError

        if self._theta is None:
            self.compute_theta()
        if self._mu is None:
            self.compute_mu()

        # Defines starting point
        if n == 0:
            x1, y1 = 0, 1
        else:
            x1 = self._points[n-1, n, 0]
            y1 = self._points[n-1, n, 1]

        mu = self._mu[n, n+1]
        theta = self._theta[n, n+1]
        k1 = brentq(lambda x: angle_between_two_lines(x, tan(theta))-mu,
                   -1E10, 0)
        x, y = intersection_of_two_rays(x1, y1, k1, 0, 0, 0)
        self._points[n, n, 0] = x
        self._points[n, n, 1] = y

    def _compute_internal_points(self, n):
        for j in xrange(n+1, self._n):
            if n == 0:
                x1, y1 = 0, 1
            else:
                x1, y1 = self._points[n-1, j]

            mu1 = self._mu[n, j+1]
            theta1 = self._theta[n, j+1]
            k1 = brentq(lambda x: angle_between_two_lines(x, tan(theta1))-mu1,
                        -1E10, 0)
            x2 = self._points[n, j-1, 0]
            y2 = self._points[n, j-1, 1]
            mu2 = self._mu[n, j]
            theta2 = self._theta[n, j]
            k2 = brentq(lambda x: angle_between_two_lines(x, tan(theta2))-mu2,
                        0, 1E10)
            x, y = intersection_of_two_rays(x1, y1, k1, x2, y2, k2)
            self._points[n, j, 0] = x
            self._points[n, j, 1] = y

    def _compute_wall_points(self, n):
        if n == 0:
            x1, y1 = 0, 1
        else:
            x1 = self._points[n-1, -1, 0]
            y1 = self._points[n-1, -1, 1]
        mu1 = self._mu[n, -1]
        theta1 = self._theta[n, -1]
        k1 = tan(theta1)

        x2 = self._points[n, -2, 0]
        y2 = self._points[n, -2, 1]
        mu2 = self._mu[n, -2]
        theta2 = self._theta[n, -2]
        k2 = brentq(lambda x: angle_between_two_lines(x, tan(theta2))-mu2,
                    0, 1E10)
        x, y = intersection_of_two_rays(x1, y1, k1, x2, y2, k2)
        self._points[n, -1, 0] = x
        self._points[n, -1, 1] = y

    def solve(self):
        if self._points is None:
            self._points = np.zeros((self._n, self.size, 2))
            self._points[:] = np.nan

        for i in xrange(self._n):
            self._compute_axis_point(i)
            self._compute_internal_points(i)
            self._compute_wall_points(i)

        return self._points

    def get_wall_points(self):
        if self._points is None:
            self.solve()

        return np.vstack((np.array([[0, 1]]),
                          self._points[:,self._n]))

    def save_plot(self, filename, dpi=None):
        if self._points is None:
            self.solve()

        fig = plt.figure()

        ax = fig.add_subplot(111, aspect='equal')
        ax.plot([0, 0, self._points[-1, -2, 0]], [1, 0, 0], 'b')

        wall_points = self.get_wall_points()
        x = wall_points[:, 0]
        y = wall_points[:, 1]
        ax.plot(x, y, 'b')

        for i in xrange(self._n):
            row = self._points[i]
            xs = row[:, 0]
            xs = xs[np.logical_not(np.isnan(xs))]
            ys = row[:, 1]
            ys = ys[np.logical_not(np.isnan(ys))]
            ax.plot(xs, ys, 'b')

        xs = self._points[:, :, 0].T
        ys = self._points[:, :, 1].T
        for i in xrange(len(xs)):
            x = xs[i]
            x = np.insert(x, 0, 0)
            x = x[np.logical_not(np.isnan(x))]
            y = ys[i]
            y = np.insert(y, 0, 1)
            y = y[np.logical_not(np.isnan(y))]
            ax.plot(x, y, 'b')

        fig.savefig(filename, dpi=dpi)
