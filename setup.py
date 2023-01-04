#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: voelzkey
"""

from setuptools import setup

setup(name='MPMQSM_preproc',
      version='0.001',
      description='preprocessing of MPMP/QSM data',
      author='YV',
      author_email='yannik.voelzke@dzne.de',
      license='MIT',
      packages=['MPMQSM_preproc'],
      entry_points={'console_scripts': ['mpmqsm_preproc   = MPMQSM_preproc.preproc:main'
                                        ]})
