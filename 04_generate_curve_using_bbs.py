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
import sys
import gmpy2

def main():

    # Test local versions of libraries

    utils.test_python_version()
    utils.test_gmpy2_version()
    utils.test_pari_version()
    utils.test_pari_seadata()
    
    now = datetime.now()
    
    # Parse command line arguments

    parser = argparse.ArgumentParser(description="Generate an Edwards curve over a given prime field, suited for cryptographic purposes.")
    parser.add_argument("input_file",
                        help="""JSON file containing the BBS parameters and the prime of the underlying field (typically, the output of
                        03_generate_prime_field_using_bbs.py.
                        """)
    parser.add_argument("output_file", help="Output file where this script will write the parameter d of the curve and the current BBS parameters.")
    parser.add_argument("--start",
                        type=int,
                        help="Number of the candidate to start with (default is 1).",
                        default=1)
    parser.add_argument("--max_nbr_of_tests",
                        type=int,
                        help="Number of candidates to test before stopping the script (default is to continue until success).")
    parser.add_argument("--fast",
                        help=""" While computing a the curve cardinality with SAE, early exit when the cardinality will obviously be divisible by
                        a small integer > 4. This reduces the time required to find the final curve, but the
                        cardinalities of previous candidates are not fully computed.
                        """,
                        default=False,
                        action="store_true")

    args = parser.parse_args()

    
    # Check arguments

    print("Checking inputs...")
    
    output_file = args.output_file
    if os.path.exists(output_file):
        utils.exit_error("The output file '%s' already exists. Exiting."%(output_file))

    input_file = args.input_file
    with open(input_file, "r") as f:
        data = json.load(f)

        
    # Declare a few important variables
        
    bbs_p = int(data["bbs_p"])
    bbs_q = int(data["bbs_q"])
    bbs_n = bbs_p * bbs_q
    bbs_s = int(data["bbs_s"]) % bbs_n
    p = int(data["p"])

    start = max(int(args.start),1)

    max_nbr_of_tests = None
    if args.max_nbr_of_tests:
        max_nbr_of_tests = int(args.max_nbr_of_tests)
        
    if not subroutines.is_strong_strong_prime(bbs_p):
        utils.exit_error("bbs_p is not a strong strong prime.")
    if not subroutines.is_strong_strong_prime(bbs_q):
        utils.exit_error("bbs_q is not a strong strong prime.")
    if not (subroutines.deterministic_is_pseudo_prime(p) and p%4 == 3):
        utils.exit_error("p is not a prime congruent to 3 modulo 4.")

        
    # Initialize BBS

    print("Initializing BBS...")
    bbs = bbsengine.BBS(bbs_p, bbs_q, bbs_s)

    
    # Info about the prime field
    
    utils.colprint("Prime of the underlying prime field:", "%d (size: %d)"%(p, gmpy2.bit_length(p)))    
    size = gmpy2.bit_length(p) # total number of bits queried to bbs for each test

    
    # Skip the first "start" candidates
    
    candidate_nbr = start-1
    bbs.skipbits(size * (start-1))


    # Start looking for "d"
    
    while True:
        
        if max_nbr_of_tests and candidate_nbr >= start + max_nbr_of_tests - 1:
            print("Did not find an adequate parameter, starting at candidate %d (included), limiting to %d candidates."%(start, max_nbr_of_tests))
            utils.exit_error("Last candidate checked was number %d."%(candidate_nbr))

        candidate_nbr += 1

        bits = bbs.genbits(size)
        d = 0
        for bit in bits:
            d = (d << 1) | bit
        print("The candidate number %d is d = %d (ellapsed time: %s)"%(candidate_nbr, d, str(datetime.now()-now)))

        
        # Test 1
        
        if not utils.check(d != 0 and d < p, "d != 0 and d < p", 1):
            continue

        # Test 2
        
        if not utils.check(gmpy2.legendre(d, p) == -1, "d is not a square modulo p", 2):
            continue
        
        # Test 3
        
        if args.fast:
            cardinality = subroutines.sea_edwards(1, d, p, 4)
        else:
            cardinality = subroutines.sea_edwards(1, d, p)
        assert(cardinality % 4 == 0)
        q = cardinality>>2
        if not utils.check(subroutines.deterministic_is_pseudo_prime(q), "The curve cardinality / 4 is prime", 3):
            continue

        # Test 4
        
        trace = p+1-cardinality
        cardinality_twist = p+1+trace
        assert(cardinality_twist % 4 == 0)
        q_twist = cardinality_twist>>2
        if not utils.check(subroutines.deterministic_is_pseudo_prime(q_twist), "The twist cardinality / 4 is prime", 4):
            continue
        
        # Test 5

        if not utils.check(q != p and q_twist != p, "Curve and twist are safe against additive transfer", 5):
            continue
        
        # Test 6

        embedding_degree = subroutines.embedding_degree(p, q)
        if not utils.check(embedding_degree > (q-1) // 100, "Curve is safe against multiplicative transfer", 6):
            continue

        # Test 7

        embedding_degree_twist = subroutines.embedding_degree(p, q_twist)
        if not utils.check(embedding_degree_twist > (q_twist-1) // 100, "Twist is safe against multiplicative transfer", 7):
            continue

        # Test 8

        D = subroutines.cm_field_discriminant(p, trace)
        if not utils.check(abs(D) >= 2**100, "Absolute value of the discriminant is larger than 2^100", 8):
            continue

        break

    
    # Find a base point

    while True:
    
        bits = bbs.genbits(size)
        y = 0
        for bit in bits:
            y = (y<<1) | bit
        u = int((1 - y**2) * gmpy2.invert(1 - d*y**2, p)) % p
        if gmpy2.legendre(u, p) == -1:
            continue
        x = gmpy2.powmod(u, (p+1) // 4, p)
        (x,y) = subroutines.add_on_edwards(x, y, x, y, d, p)
        (x,y) = subroutines.add_on_edwards(x, y, x, y, d, p)
        if (x, y) == (0, 1):
            continue

        assert((x**2 + y**2) % p == (1 + d*x**2*y**2) % p)
        
        break

    
    # Print some informations
    
    utils.colprint("Number of the successful candidate:", str(candidate_nbr))
    utils.colprint("Edwards elliptic curve parameter d is:", str(d))
    utils.colprint("Number of points:", str(cardinality))
    utils.colprint("Number of points on the twist:", str(cardinality_twist))
    utils.colprint("Embedding degree of the curve:", "%d"%embedding_degree)
    utils.colprint("Embedding degree of the twist:", "%d"%embedding_degree_twist)
    utils.colprint("Discriminant:", "%d"%D)
    utils.colprint("Trace:", "%d"%trace)
    utils.colprint("Base point coordinates:", "(%d, %d)"%(x, y))

    
    # Save p, d, x, y, etc. to the output_file

    print("Saving the parameters to %s"%output_file)
    bbs_s = bbs.s
    with open(output_file, "w") as f:
        json.dump({"p": int(p),
                   "bbs_p": int(bbs_p),
                   "bbs_q": int(bbs_q),
                   "bbs_s": int(bbs_s),
                   "candidate_nbr": int(candidate_nbr),
                   "d": int(d),
                   "cardinality": cardinality,
                   "cardinality_twist": cardinality_twist,
                   "embedding_degree": embedding_degree,
                   "embedding_degree_twist": embedding_degree_twist,
                   "discriminant": D,
                   "trace": trace,
                   "base_point_x": x,
                   "base_point_y": y},
                  f,
                  sort_keys=True)

    
if __name__ == "__main__":
    main()
