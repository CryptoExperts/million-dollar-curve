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

import sys
import platform
import os.path
import subroutines
import gmpy2


def exit_error(s):
    sys.exit("[ERROR] %s"%s)

    
def colprint(col1,col2 = "",col1_size = 50):
    col1_size = max(col1_size, len(col1))
    if col2 != "":
        print(col1 + " " * (col1_size - len(col1)) + col2)
    else:
        print(col1)

        
def check(test, test_description="", test_number=None):
    
    if test_description != "":
        if test_number != None:
            print("\tTest %d: %s... "%(test_number,test_description), end="", flush=True)
        else:
            print("\tTesting %s... "%(test_description), end="", flush=True)
            
    if test:
        print("Success")
        return True
    else:
        print("Failure")
        return False
    
    
def test_gmpy2_version():
    expected_gmpy2_version = '2.0.7'
    local_gmpy2_version = gmpy2.version()
    if local_gmpy2_version != expected_gmpy2_version:
        print('[WARNING] You are using version %s of gmpy2. These scripts have been tested with version %s'
              %(local_gmpy2_version, expected_gmpy2_version), file=sys.stderr)

        
def test_python_version():
    expected_python_version = '3.4.3'
    local_python_version = platform.python_version()
    if local_python_version != expected_python_version:
        print('[WARNING] You are using version %s of Python. These scripts have been tested with version %s'
              %(local_python_version, expected_python_version), file=sys.stderr)

        
def test_pari_version():
    expected_pari_version = '2.7.4'
    local_pari_version = subroutines.pari_version()
    if local_pari_version != expected_pari_version:
        print('[WARNING] You are using version %s of PARI. These scripts have been tested with version %s'
              %(local_pari_version, expected_pari_version), file=sys.stderr)
    
        
def test_pari_seadata():
    """Look for the seadata package for PARI. Exit if it cannot be found."""
    datadir = subroutines.pari_cfg_datadir()
    seadata_path = os.path.join(datadir,"seadata")
    if not os.path.exists(seadata_path):
        exit_error("Cannot find the seadata.tgz package for PARI. Please install this package in %s. See http://pari.math.u-bordeaux.fr/packages.html."%(datadir))
