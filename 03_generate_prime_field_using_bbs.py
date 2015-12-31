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
import subroutines
import utils
import gmpy2

def main():

    # Test local versions of libraries

    utils.test_python_version()
    utils.test_gmpy2_version()
    
    # Parse command line arguments

    parser = argparse.ArgumentParser(description="Generate a prime field, suited for being the underlying field of a twist-secure Edwards curve.")
    parser.add_argument("input_file", help="JSON file containing the BBS parameters (typically, the output of 02_generate_bbs_parameters.py).")
    parser.add_argument("output_file", help="Output file where this script will write the prime of the field and the current BBS parameters.")
    parser.add_argument("prime_size", type=int, help="Size of the prime (e.g. 256 bits)")
    
    args = parser.parse_args()

    
    # Check arguments

    output_file = args.output_file
    if os.path.exists(output_file):
        utils.exit_error("The output file '%s' already exists. Exiting."%(output_file))

    size = int(args.prime_size)

    input_file = args.input_file
    with open(input_file, "r") as f:
        data = json.load(f)        
    bbs_p = int(data["bbs_p"])
    bbs_q = int(data["bbs_q"])
    bbs_n = bbs_p * bbs_q
    bbs_s = int(data["bbs_s"]) % bbs_n

    
    # Check inputs

    print("Checking inputs...")
    if not subroutines.is_strong_strong_prime(bbs_p):
        utils.exit_error("bbs_p is not a strong strong prime.")
    if not subroutines.is_strong_strong_prime(bbs_q):
        utils.exit_error("bbs_q is not a strong strong prime.")

        
    # Initialize BBS

    bbs = bbsengine.BBS(bbs_p, bbs_q, bbs_s)

    
    # generate a "size"-bit prime "p"

    candidate_nbr = 0
    print("Generating a prime field Fp (where p is congruent to 3 mod 4)...")
    while True:
        candidate_nbr += 1
        bits = [1] + bbs.genbits(size-3) + [1,1]
        assert(len(bits) == size)
        p = 0
        for bit in bits:
            p = (p << 1) | bit
        assert(p % 4 == 3)
        assert(gmpy2.bit_length(p) == size)
        if subroutines.deterministic_is_pseudo_prime(p):
            break
    utils.colprint("%d-bit prime found:"%size, str(p))
    utils.colprint("The good candidate was number: ", str(candidate_nbr))

    
    # Save p and the current bbs parameters to the output_file

    print("Saving p and the BBS parameters to %s"%(output_file))
    bbs_s = bbs.s
    with open(output_file, "w") as f:
        json.dump({"p": int(p), 
                   "bbs_p": int(bbs_p), 
                   "bbs_q": int(bbs_q), 
                   "bbs_s": int(bbs_s)}, 
                  f,
                  sort_keys=True)

    
if __name__ == "__main__":
    main()
