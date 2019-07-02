#!/usr/bin/env python

# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'ifoche@cnb.csic.es'
# *
# **************************************************************************

import os
from os.path import join
from glob import glob
from datetime import datetime


Import('env')


AddOption('--no-opencv', dest='opencv', action='store_false', default=True,
          help='Avoid compilation of opencv programs')
AddOption('--no-scipy', dest='scipy', action='store_false', default=True,
          help='Avoid compilation with scipy support')



# Define some variables used by Scons. Note that some of
# the variables will be passed by Scipion in the environment (env).

env['CUDA_SDK_PATH'] = os.environ.get('CUDA_SDK_PATH', '')
env['CUDA_LIB_PATH'] = os.environ.get('CUDA_LIB_PATH', '')

get = lambda x: os.environ.get(x, '0').lower() in ['true', 'yes', 'y', '1']

gtest = get('GTEST')
debug = get('DEBUG')

# Read some flags
CYGWIN = env['PLATFORM'] == 'cygwin'
MACOSX = env['PLATFORM'] == 'darwin'
MINGW = env['PLATFORM'] == 'win32'

XMIPP_PATH = Dir('.').abspath
XMIPP_BUNDLE = Dir('..').abspath


#  ***********************************************************************
#  *                      Xmipp C++ Libraries                            *
#  ***********************************************************************

# Create a shortcut and customized function
# to add the Xmipp CPP libraries
def addLib(name, **kwargs):
    # Install all libraries in scipion/software/lib
    # COSS kwargs['installDir'] = '#software/lib'
    # Add always the xmipp path as -I for include and also xmipp/libraries
    incs = kwargs.get('incs', [])
    kwargs['incs'] = incs

    deps = kwargs.get('deps', [])
    kwargs['deps'] = deps

    # Add libraries in libs as deps if not present
    libs = kwargs.get('libs', [])
    for lib in libs:
        deps.append(lib)

    # If pattern not provided use *.cpp as default
    patterns = kwargs.get('patterns', '*.cpp')
    kwargs['patterns'] = patterns
    lib = env.AddCppLibrary(name, **kwargs)
    	
    env.Alias('xmipp-libs', lib)

    return lib


# Gtest
#addLib('XmippGtest',
#       dirs=['external'],
#       patterns=['gtest/*.cc'],
#       default=gtest,
#       libs=['pthread']
#       )

def getHdf5Name(libdirs):
    for dir in libdirs:
        if os.path.exists(os.path.join(dir.strip(),"libhdf5.so")):
            return "hdf5"
        elif os.path.exists(os.path.join(dir.strip(),"libhdf5_serial.so")):
            return "hdf5_serial"
    return "hdf5"

# Data
addLib('XmippCore',
       patterns=['*.cpp','*.c','bilib/*.cc','alglib/*.cpp', 'utils/*.cpp'],
       dirs=['core'] * 5, # one relative path for each pattern
       libs=['fftw3', 'fftw3_threads',
             getHdf5Name(env['EXTERNAL_LIBDIRS']),'hdf5_cpp',
             'tiff',
             'jpeg',
             'sqlite3',
             'pthread'])

# Python binding
def remove_prefix(text, prefix):
    return text[text.startswith(prefix) and len(prefix):]
env['PYTHONINCFLAGS'] = os.environ.get('PYTHONINCFLAGS', '').split()
if len(env["PYTHONINCFLAGS"])>0:
    python_incdirs = [remove_prefix(os.path.expandvars(x),"-I") for x in env["PYTHONINCFLAGS"]]
else:
    python_incdirs = []

addLib('xmippCore.so',
       dirs=['bindings'],
       patterns=['python/*.cpp'],
       incs=python_incdirs,
       libs=['python2.7', 'XmippCore'],
       prefix='', target='xmippCore')


#  ***********************************************************************
#  *                      Xmipp Scripts                                  *
#  ***********************************************************************

XmippAlias = env.Alias('xmipp', ['xmipp-libs'])
Return('XmippAlias')
