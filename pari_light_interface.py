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

import os
import ctypes
import ctypes.util

libpari = ctypes.CDLL(ctypes.util.find_library("pari"))


_fx_pari_version          = libpari.pari_version
_fx_pari_version.argtypes = None
_fx_pari_version.restype  = ctypes.c_void_p

def pari_version():
    return _fx_pari_version()


_fx_pari_sd_datadir          = libpari.sd_datadir
_fx_pari_sd_datadir.argtypes = [ctypes.c_char_p, ctypes.c_long]
_fx_pari_sd_datadir.restype  = ctypes.c_void_p

def pari_sd_datadir():
    return _fx_pari_sd_datadir(None, 3) # Deduced from paridecl.h


_fx_pari_init          = libpari.pari_init
_fx_pari_init.argtypes = [ctypes.c_size_t, ctypes.c_ulong]
_fx_pari_init.restype  = None

def pari_init(size, maxprime):
    _fx_pari_init(size, maxprime)


_fx_pari_close          = libpari.pari_close
_fx_pari_close.argtypes = None
_fx_pari_close.restype  = None

def pari_close():
    _fx_pari_close()
    

_fx_pari_gp_read_str          = libpari.gp_read_str
_fx_pari_gp_read_str.argtypes = [ctypes.c_char_p]
_fx_pari_gp_read_str.restype  = ctypes.c_void_p

def pari_gp_read_str(s):
    return _fx_pari_gp_read_str(str(s).encode("UTF-8"))


_fx_pari_GENtostr          = libpari.GENtostr
_fx_pari_GENtostr.argtypes = [ctypes.c_void_p]
_fx_pari_GENtostr.restype  = ctypes.c_char_p

def pari_GENtostr(a):
    return _fx_pari_GENtostr(a)


_fx_pari_addii          = libpari.addii
_fx_pari_addii.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
_fx_pari_addii.restype  = ctypes.c_void_p

def pari_addii(a, b):
    return _fx_pari_addii(a, b)


_fx_pari_Fp_ellcard_SEA          = libpari.Fp_ellcard_SEA
_fx_pari_Fp_ellcard_SEA.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_long]
_fx_pari_Fp_ellcard_SEA.restype  = ctypes.c_void_p

def pari_Fp_ellcard_SEA(a4, a6, p, s):
    return _fx_pari_Fp_ellcard_SEA(a4, a6, p, s)


_fx_pari_Z_factor          = libpari.Z_factor
_fx_pari_Z_factor.argtypes = [ctypes.c_void_p]
_fx_pari_Z_factor.restype  = ctypes.c_void_p

def pari_Z_factor(n):
    return _fx_pari_Z_factor(n)


def pari_gel(x, i):
    s = ctypes.sizeof(ctypes.c_void_p)
    v = ctypes.c_void_p.from_address(x + i*s)
    return v.value


def pari_lg(z):
    s = ctypes.sizeof(ctypes.c_void_p)
    header = ctypes.c_long.from_address(z).value
    mask = (1 << (8*s - 8)) - 1 # Deduced from parigen.h
    return header & mask
