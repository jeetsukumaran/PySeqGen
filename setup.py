#! /usr/bin/env python

############################################################################
## setup.py
##
##  Part of the PySeqGen wrapper library.
##
##  Copyright 2007 Jeet Sukumaran and Mark T. Holder.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
############################################################################

import ez_setup
import sys
ez_setup.use_setuptools()

from setuptools import setup, find_packages, Extension

setup(name = "pyseqgen",
      version = "2.00",
      packages = find_packages(),
      author="Jeet Sukumaran and Mark T. Holder",
      author_email='jeetsukumaran@gmail.com',
      maintainer = "Jeet Sukumaran and Mark T. Holder",
      maintainer_email = "jeetsukumaran@gmail.com mholder@ku.edu",
      description = "Wrapper around Seq-Gen",
      license='LGPL 3+',
      install_requires=[
          "DendroPy >= 2.5.0",
      ],
      long_description="""\
Simple wrappers around Seq-Gen molecular sequence generator""",
      classifiers = [
            "Development Status :: 2 - Pre-Alpha",
            "Environment :: Console",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU Library or  General Public License (GPL)",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            ]
    )
