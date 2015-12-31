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

import gmpy2

class BBS:
    """Blum Blum Shub Pseudo Random Generator Engine"""

    def __init__(self, p, q, s, shift=0):
        self.p = p
        self.q = q            
        self.n = p * q
        power = (2**shift) % ((p-1) * (q-1))
        self.s = s % self.n
        self.s = gmpy2.powmod(self.s, power, self.n)

    def genbit(self):
        self.s = (self.s**2) % self.n
        return self.s & 1

    def genbits(self, k):
        bits = []
        for i in range(k):
            bits.append(self.genbit())
        return bits

    def skipbits(self,k):
        power = (2**k) % ((self.p-1) * (self.q-1))
        self.s = gmpy2.powmod(self.s, power, self.n)
