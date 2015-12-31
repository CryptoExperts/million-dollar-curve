#!/usr/bin/env python3

# This file is part of Million Dollar Curve

# Copyright (C) 2015, 2016  CryptoExperts

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

import pari_light_interface
import gmpy2

FIRST_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
                103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211,
                223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337,
                347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461,
                463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
                607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739,
                743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881,
                883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021,
                1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129,
                1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277,
                1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409,
                1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511,
                1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619]

def pari_version():

    pari_light_interface.pari_init(100000000, 0)

    _v = pari_light_interface.pari_version()
    _x = pari_light_interface.pari_gel(_v, 1)
    _y = pari_light_interface.pari_gel(_v, 2)
    _z = pari_light_interface.pari_gel(_v, 3)
    x = int(pari_light_interface.pari_GENtostr(_x))
    y = int(pari_light_interface.pari_GENtostr(_y))
    z = int(pari_light_interface.pari_GENtostr(_z))

    pari_light_interface.pari_close()
    
    return str("%d.%d.%d"%(x,y,z))

def pari_cfg_datadir():

    pari_light_interface.pari_init(100000000, 0)

    _s = pari_light_interface.pari_sd_datadir()
    s = str(pari_light_interface.pari_GENtostr(_s).decode('utf-8'))[1:-1]
    
    pari_light_interface.pari_close()

    return s
    
def sea_weierstrass(a, b, p, s=0):

    pari_light_interface.pari_init(100000000, 0)

    _a = pari_light_interface.pari_gp_read_str(str(a))
    _b = pari_light_interface.pari_gp_read_str(str(b))
    _p = pari_light_interface.pari_gp_read_str(str(p))

    _t = pari_light_interface.pari_Fp_ellcard_SEA(_a, _b, _p, s)
    t  = pari_light_interface.pari_GENtostr(_t)
    
    pari_light_interface.pari_close()
    
    return int(t)

def sea_montgomery(A, B, p, s=0):
    (a,b) = _weierstrass_parameters_from_montgomery_parameters(A, B, p)
    return sea_weierstrass(a, b, p, s)

def sea_edwards(a, d, p, s=0):
    (A,B) = _montgomery_parameters_from_edwards_parameters(a, d, p)
    return sea_montgomery(A, B, p, s)

def _weierstrass_parameters_from_montgomery_parameters(A, B, p):
    a = ((3 - A**2) * gmpy2.invert(3*B**2, p)) % p
    b = ((2*A**3 - 9*A) * gmpy2.invert(27*B**3, p)) % p
    return (a, b)

def _montgomery_parameters_from_edwards_parameters(a, d, p):
    A = (2 * (a+d) * gmpy2.invert(a - d, p)) % p
    B = (4 * gmpy2.invert(a - d, p)) % p
    return (A, B)

def factor(n):

    pari_light_interface.pari_init(100000000, 0)

    _n = pari_light_interface.pari_gp_read_str(str(n))

    _A = pari_light_interface.pari_Z_factor(_n)
    _P = pari_light_interface.pari_gel(_A, 1)
    _M = pari_light_interface.pari_gel(_A, 2)
    l = pari_light_interface.pari_lg(_P)
    
    f = []
    for i in range(1, l):
        _p = int(pari_light_interface.pari_gel(_P, i))
        p = int(pari_light_interface.pari_GENtostr(_p))
        _m = int(pari_light_interface.pari_gel(_M, i))
        m = int(pari_light_interface.pari_GENtostr(_m))
        if f and len(f) > 0:
            assert(p > f[-1][0]) # if this fails, add some code that makes sure f is sorted
        f.append([p, m])

    pari_light_interface.pari_close()
    
    return f

def cm_field_discriminant(p, t):

    # Compute s^2, the largest square dividing t^2-4p
    s_square = 1
    factors = factor(t**2 - 4*p)
    for f in factors:
        if f[1] % 2 == 0:
            s_square *= f[0]**f[1]
        else:
            s_square *= f[0]**(f[1]-1)
    assert((t**2 - 4*p) % s_square == 0)

    # Compute D
    D = (t**2 - 4*p) // s_square
    if D % 4 != 1:
        D *= 4

    return D

def embedding_degree(p, q):
    """Given p (typically the prime of the base field) and q (typically the order of a large subgroup the EC), return the
    embedding degree of the curve, i.e., the smallest m such that p^m = 1 (mod q).
    
    Keyword arguments:
    p -- prime of the base field
    q -- the order of a large subgroup the EC
    """
    m = q - 1
    factors = factor(m)
    for f in factors:
        for i in range(f[1]):
            if (p^(m // f[0])) % q == 1:
                m = m // f[0]

    return m


def is_strong_prime(p):
    if not deterministic_is_pseudo_prime(p):
        return False
    p = (p-1) >> 1
    return deterministic_is_pseudo_prime(p)


def is_strong_strong_prime(p):
    if not deterministic_is_pseudo_prime(p):
        return False
    p = (p-1) >> 1
    if not deterministic_is_pseudo_prime(p):
        return False
    p = (p-1) >> 1
    return deterministic_is_pseudo_prime(p)


def is_strong_strong_prime_generator(p):
    if not deterministic_is_pseudo_prime(p):
        return False
    p = (p<<1) + 1
    if not deterministic_is_pseudo_prime(p):
        return False
    p = (p<<1) + 1
    return deterministic_is_pseudo_prime(p)

def add_on_edwards(x1, y1, x2, y2, d, p):
    x = int((x1*y2 + x2*y1) * gmpy2.invert(1+d*x1*x2*y1*y2, p)) % p
    y = int((y1*y2 - x1*x2) * gmpy2.invert(1-d*x1*x2*y1*y2, p)) % p
    return (x,y)

def deterministic_is_pseudo_prime(n, k=64):
    assert(k <= len(FIRST_PRIMES))
    if n in FIRST_PRIMES:
        return True
    if n < max(FIRST_PRIMES):
        return False
    if n & 1 == 0:
        return False
    # Find s and t
    s = 0
    t = n - 1
    while t & 1 == 0:
        s += 1
        t = t >> 1
    assert(n == 2**s * t + 1)
    # main loop
    for j in range(k):
        b = FIRST_PRIMES[j]
        x = gmpy2.powmod(b, t, n)
        i = 0
        if x != 1:
            while x != n - 1:
                x = gmpy2.powmod(x, 2, n)
                i += 1
                if i == s or x == 1:
                    return False
    return True
