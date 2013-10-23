from __future__ import division

__author__ = 'psinger'

"""
Functions for calculating the statistical significant differences between two dependent or independent correlation
coefficients.
Adopted from the R package http://personality-project.org/r/html/paired.r.html
"""

import numpy as np
from scipy.stats import t, norm

def dependent_corr(xy, xz, yz, n, twotailed = True):
    """
    @param xy: correlation coefficient between x and y
    @param xz: correlation coefficient between x and z
    @param yz: correlation coefficient between y and z
    @param n: number of elements in x, y and z
    @param twotailed: whether to calculate a one or two tailed test
    @return: t and p-val
    """
    d = xy - xz
    determin = 1 - xy * xy - xz * xz - yz * yz + 2 * xy * xz * yz
    av = (xy + xz)/2
    cube = (1 - yz) * (1 - yz) * (1 - yz)

    t2 = d * np.sqrt((n - 1) * (1 + yz)/(((2 * (n - 1)/(n - 3)) * determin + av * av * cube)))
    p = 1 - t.cdf(abs(t2), n - 2)

    if twotailed:
        p *= 2

    return t2, p

def independent_corr(xy, xz, n, n2 = None, twotailed = True):
    """
    @param xy: correlation coefficient between x and y
    @param xz: correlation coefficient between x and z
    @param n: number of elements in xy
    @param n2: number of elements in xz (if distinct from n)
    @param twotailed: whether to calculate a one or two tailed test
    @return: z and p-val
    """
    xy_z = 0.5 * np.log((1 + xy)/(1 - xy))
    xz_z = 0.5 * np.log((1 + xz)/(1 - xz))
    if n2 is None:
        n2 = n

    se_diff_r = np.sqrt(1/(n - 3) + 1/(n2 - 3))
    diff = xy_z - xz_z
    z = abs(diff / se_diff_r)
    p = (1 - norm.cdf(z))
    if twotailed:
        p *= 2

    return z, p

print dependent_corr(.396, .179, .088,200)
print independent_corr(.560, .588, 100, 353, twotailed=False)

