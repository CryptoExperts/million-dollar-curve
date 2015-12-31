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
import re
import math
from datetime import datetime,timedelta
import gmpy2
import utils


def main():

    # Test local versions of libraries

    utils.test_python_version()
    utils.test_gmpy2_version()

    # Parse command line arguments
    
    parser = argparse.ArgumentParser(description="From a table of draws, output a seed of appropriate length.")

    parser.add_argument("input_draw_file", help="""Text file with draws.""")
    parser.add_argument("output_seed_file", help="""JSON file where we can store the seed computed from the draws.""")
    parser.add_argument("entropy_to_gather", help="""Minimum entropy to extract before drawing lone bits.""")
    parser.add_argument("--nbr_lone_bits", type=int, help="""Number of lone bits to extract.""", default=0)
    
    args = parser.parse_args()


    # Check arguments
    
    output_seed_file = args.output_seed_file
    if os.path.exists(output_seed_file):
        utils.exit_error("The output file '%s' already exists. Exiting."%(output_seed_file))


    # Declare a few important variables
    
    two_pow_entropy_to_gather = (1<<int(args.entropy_to_gather))
    
    seed = 0
    L = 1 # before lone bits are drawn, seed lies in [0,L - 1]

    lone_bits_part = 0
    nbr_lone_bits = args.nbr_lone_bits


    # Scan the input file, construct the seed
    
    with open(args.input_draw_file, "r") as f:
        
        for line in f:

            if not line or line.strip() == "" or line.startswith("#"):
                continue

            (draw_id, m, n, draw) = re.split("\s+", line.strip(), maxsplit=3)

            if draw == "None":
                continue
            
            m = int(m)
            n = int(n)
            draw = [ int(x) for x in draw.split(",") ]
            index = index_from_draw(draw,m)

            if L < two_pow_entropy_to_gather:
                
                print("Draw %s used to extract entropy"%(draw_id))
                seed = gmpy2.bincoef(n,m)*seed + index
                L *= gmpy2.bincoef(n,m)
                
            else:
                
                print("Draw %s used to extract a lone bit"%(draw_id))
                b = index & 1
                seed += L * (b << (args.nbr_lone_bits - nbr_lone_bits))
                nbr_lone_bits -= 1

            if L >= two_pow_entropy_to_gather and nbr_lone_bits == 0:
                break

    if nbr_lone_bits > 0 or L < two_pow_entropy_to_gather:
        utils.exit_error("There wasn't enough draws to collect to request quantity of entropy and lone bits.")

    seed_upper_bound = L * 2**(args.nbr_lone_bits)
    seed_entropy = math.floor(gmpy2.log2(seed_upper_bound))
    print("The seed contains more than %d bits of entropy (including the %s lone bits)."%(seed_entropy,args.nbr_lone_bits))
    print("The seed is %d"%(seed))

    print("Saving the seed to %s"%(output_seed_file))
    with open(output_seed_file, "w") as f:
        json.dump({"seed": int(seed),
                   "seed_upper_bound": int(seed_upper_bound),
                   "approx_seed_entropy": int(seed_entropy),
                   "lone_bits": int(args.nbr_lone_bits)}, 
                  f,
                  sort_keys=True)

            
def index_from_draw(draw,m):
    index = 0
    draw = sorted(draw)
    for i in range(m):
        index += gmpy2.bincoef(draw[i]-1, i+1)
    return index


if __name__ == "__main__":
    main()
