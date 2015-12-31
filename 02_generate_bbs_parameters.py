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
import json
import os
import subroutines
import utils
import math
import gmpy2

def main():

    # Test local versions of libraries

    utils.test_python_version()
    utils.test_gmpy2_version()

    # Parse command line arguments
    
    parser = argparse.ArgumentParser(description="Generate BBS parameters.")
    
    parser.add_argument("input_file", help="""JSON file containing the seed used for generating the pseudo strong 
                                              strong prime (the name is "seed"). The required
                                              quantity of entropy it should contain depends on bitsize. As a rule of
                                              thumb the seed should contain at least 4*bitsize bits of entropy.""")
    parser.add_argument("output_file", help="""Output JSON file where this script will write the two generated strong
                                               strong primes "p" and "q". The output file should not exist already.""")
    parser.add_argument("min_prime_bitsize", type=int, help="minimum strong strong prime bit size (e.g. 2048).")
    
    args = parser.parse_args()

    
    # Check arguments
    
    output_file = args.output_file
    if os.path.exists(output_file):
        utils.exit_error("The output file '%s' already exists. Exiting."%(output_file))


    # Declare a few important variables
        
    min_prime_bitsize = args.min_prime_bitsize

    input_file = args.input_file
    with open(input_file, "r") as f:
        data = json.load(f)        
    seed = int(data["seed"])
    seed_upper_bound = int(data["seed_upper_bound"])
    approx_seed_entropy = math.floor(gmpy2.log2(seed_upper_bound))

    utils.colprint("Minimum strong strong prime size:", str(min_prime_bitsize))
    utils.colprint("Approximate seed entropy:", str(approx_seed_entropy))

    
    # Precomputations

    first_primes = [2]                     # List of the first primes
    PI = 2                                 # Product of the primes in "first_primes"
    strong_strong_integers = [[1]]         # strong_strong_integers[i] is the list of all strong strong integers modulo
                                           # first_primes[i]
    number_of_strong_strong_integers = [1] # number_of_strong_strong_integers[i] is the number of elements of the list
                                           # strong_strong_integers[i]
    C = 1                                  # Product of the elements of "number_of_strong_strong_integers"
    
    while not 2**(min_prime_bitsize-2) < PI:
        p = int(gmpy2.next_prime(first_primes[-1]))
        first_primes.append(p)
        PI *= p
        ssi = [c for c in range(p) if is_strong_strong_basis(c, p)]
        strong_strong_integers.append(ssi)
        number_of_strong_strong_integers.append(len(ssi))
        C *= len(ssi)

    utils.colprint("Number of primes considered:", str(len(first_primes)))
    utils.colprint("Number of strong strong integers to choose from:", "about 2^%f"%(gmpy2.log2(C)))

    
    # Check that the seed is long enough

    if seed_upper_bound < C**2 * (1 << (2 * min_prime_bitsize)):
        utils.exit_error("The seed does not contain the required entropy.")

        
    # Precomputations for the CRT

    mu    = [gmpy2.divexact(PI,p) for p in first_primes]
    delta = [gmpy2.invert(x,y) for x,y in zip(mu,first_primes)]
    gamma = [gmpy2.mul(x,y) for x,y in zip(mu,delta)]


    # Generate the first strong prime
    
    print("Generating the first strong strong prime...")
    (p,seed) = generate_strong_strong_prime(seed,
                                            min_prime_bitsize,
                                            strong_strong_integers,
                                            number_of_strong_strong_integers,
                                            gamma,
                                            PI)
    utils.colprint("\tThis is the first strong strong prime:", str(p))

    
    # Generate the second strong prime
    
    print("Generating the second strong strong prime...")
    (q,seed) = generate_strong_strong_prime(seed,
                                            min_prime_bitsize,
                                            strong_strong_integers,
                                            number_of_strong_strong_integers,
                                            gamma,
                                            PI)
    utils.colprint("\tThis is the second strong strong prime:", str(q))

    
    # Generate the BBS start

    print("Generating the BBS starting point...")    
    n = p*q
    s = seed % n
    while s == 0 or s == 1 or s == p or s == q:
        s = (s+1) % n
    s0 = (s**2) % n
    utils.colprint("\tThis is the starting point s0 of BBS:", str(s0))

    
    # Save p,q, and s to the output_file

    print("Saving p,q, and s0 to %s"%(output_file))
    with open(output_file, "w") as f:
        json.dump({"bbs_p": int(p), 
                   "bbs_q": int(q), 
                   "bbs_s": int(s0)}, 
                  f,
                  sort_keys=True)


    
def generate_strong_strong_prime(seed, min_bitsize,strong_strong_integers,number_of_strong_strong_integers,gamma,PI):
    """Return a strong strong prime deterministically determined from the input parameters, and what remains of the seed.

    Depending on the target prime size "min_bitsize", we need to find the appropriate table of first primes
    [p_0,...,p_{f-2},p_{f-1}] such that PI = p_0 * ... * p_{f-1} is larger than 2**(min_bitsize-2), but such that p_0 *
    ... * p_{f-2} is not. We'll store the primes in a table called "first_primes", such that first_primes[i] = p_i. The
    variable "strong_strong_integers" will be a list of lists, such that strong_strong_integers[i] is the list of
    integers r such that r, 2r+1 and 2(2r+1)+1 are invertible modulo first_primes[i]. Taking the inverse CRT on any
    
    [c_0,c_1,...,c_{f-1}] = [strong_strong_integers[0][i],strong_strong_integers[1][j],...,strong_strong_integers[f-1][k]]

    gives an integer c such that c, 2c+1, and 4c+3 are invertible modulo PI, which makes c a good candidate for being
    a strong strong prime generator (i.e., 4c+3 a good candidate for being a strong strong prime). The seed allows to
    determine in initial array [c_0,c_1,...,c_{f-1}] and thus an initial candidate c. Going from one such array to the
    other is done deterministically.

    """

    # Consume the seed and update it
    
    (indexes, seed) = list_of_indexes_from_seed(seed, number_of_strong_strong_integers)
    
    candidate_nbr = 0
    while True:

        alpha = [x[i] for x, i in zip(strong_strong_integers, indexes)] # alpha is in Z2* x Z3* x Z5* x ..... Apply the inverse CRT
        c = sum([x*y for x, y in zip(alpha, gamma)]) % PI
        candidate_nbr += 1
        
        if gmpy2.bit_length(c) >= min_bitsize-2 and subroutines.is_strong_strong_prime_generator(c):
            break

        indexes = next_indexes(indexes, number_of_strong_strong_integers)

    print("\tThe successful candidate is the number %d"%(candidate_nbr))

    return (4*c + 3,seed)

    

def is_strong_strong_basis(alpha, p):
    """Return True if alpha, 2*alpha+1, and 2*(2*alpha+1) + 1 are invertible modulo p,
    and false otherwise.

    """

    if (alpha % p) == 0:
        return False
    if ((2*alpha + 1) % p) == 0:
        return False
    if ((2*(2*alpha + 1) + 1) % p) == 0:
        return False
    return True


def list_of_indexes_from_seed(seed, max_indexes):
    """Return a list of len(max_indexes) indexes computed from the seed, such that
    0 <= indexes[i] < max_indexes[i]. Also return what remains from the seed.
    """
    indexes = []
    for i in range(len(max_indexes)):
        r = seed % max_indexes[i]
        seed = (seed-r) // max_indexes[i]
        indexes.append(r)
    return (indexes,seed)


def next_indexes(indexes, max_indexes):
    i = 0
    while True:
        indexes[i] = (indexes[i]+1) % max_indexes[i]
        if indexes[i] == 0:
            i = (i+1) % len(indexes)
        else:
            break
    return indexes

    
if __name__ == "__main__":
    main()
