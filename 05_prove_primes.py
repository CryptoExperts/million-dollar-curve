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

import argparse
import bbsengine
import json
import os
import utils
import subroutines
from datetime import datetime
import math
import gmpy2
import sys

def main():

    now = datetime.now()
    
    # Parse command line arguments

    parser = argparse.ArgumentParser(description=
                                     """The script takes as an input a list of integers. If any of those integers fails 
                                     a pseudo primality test, this script exists immediately. Otherwise, it tries to 
                                     output a proof of primality for each of these pseudo primes.
                                     """)
    
    parser.add_argument("integers", type=int, nargs="+", help="List of all integers to consider.")

    args = parser.parse_args()


    # Check arguments
    
    for n in args.integers:
        if not subroutines.deterministic_is_pseudo_prime(n):
            utils.exit_error("%d is not prime."%(n))
        
    # Declare a few important variables. In particular, large_factors[p] will contain a list [[p1,m1],[p2,m2],...]  such
    # that p1^m1*p2^m2*... > sqrt(p), for all "p" in "pseudo_primes".

    pseudo_primes = set(args.integers)
    large_factors = {} 

    
    # Compute the dictionnary "large_factors"

    while True:
        
        if not pseudo_primes:
            break
        
        p = max(pseudo_primes)
        pseudo_primes.remove(p)

        if p == 2:
            continue
        
        print("Factoring %d - 1"%(p))
        
        all_factors = subroutines.factor(p-1)

        A = 1
        f = []
        while A <= math.sqrt(p):
            [q,m] = all_factors.pop()
            A *= q**m
            f += [[q,m]]
            pseudo_primes.add(q)
            
        large_factors[p] = list(reversed(f))

    # Prove primes
        
    proven_primes = {2: []} # For N > 2, proven_primes[N] will be an array [large_factors[N],a], where proof is a
                            # dictionnary s.t. len(a) == len(large_factors[N]) and a[p] is the a_p corresponding to the
                            # factor p = large_factors[N][p] in the Pocklington method.

    while True:

        if not large_factors:
            break

        N = min(large_factors.keys())
        # Generalized Pocklington method to show that N is prime
        f = large_factors.pop(N) # large factors of N - 1
        a = {}
        for p,m in f:
            assert((N-1) % p == 0)
            for a_p in range(2,N):
                if gmpy2.powmod(a_p, N-1, N) != 1:
                    continue
                if gmpy2.gcd(gmpy2.powmod(a_p, (N-1)//p, N) - 1, N) != 1:
                    continue
                break
            a[p] = a_p
        proven_primes[N] = [f,a]

        
    # Print proofs

    for N in sorted(proven_primes.keys()):
        if N == 2:
            continue
        print("Proof that N = %d is prime:"%(N))
        f = proven_primes[N][0]
        a = proven_primes[N][1]
        A = 1
        for p,m in proven_primes[N][0]:
            A *= p**m
        assert((N-1) % A == 0)
        B = (N-1) // A
        print("\tN - 1 = A * B with")
        print("\tA = %d = %s"%(A,factors_to_string(f)))
        print("\tB = %d"%(B))
        assert(gmpy2.gcd(A,B) == 1)
        assert(A > math.sqrt(N))
        print("\tA and B are relatively prime and A > sqrt(N).")
        print("\tPrime factor(s) of A: %s"%(", ".join([str(p) for p,m in f])))
        for p,m in f:
            assert(gmpy2.powmod(a[p], N-1, N) == 1)
            assert(gmpy2.gcd(gmpy2.powmod(a[p], (N-1) // p, N), N) == 1)
            print("\tFor p = %d, we have %d^(N-1) mod N = 1 and gcd(%d^((N-1)/p) - 1, N) = 1"%(p, a[p], a[p]))


def factors_to_string(f):
    s = ""
    for p,m in f:
        s += "%d^%d * "%(p,m)
    return s[0:-3]


if __name__ == "__main__":
    main()
