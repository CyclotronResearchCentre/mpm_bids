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
      package_data={'MPMQSM_preproc': ['EP3D_mtsaturation.ini', 'EP3D_pushmt_yv.ini', 'spm12/**/*', 'hMRI-toolbox/**/*', 'call_batch.m', 'hmri_defaults_SS_scaifield.m']},
      include_package_data=True,
      entry_points={'console_scripts': ['mpmqsm_preproc   = MPMQSM_preproc.preproc:main'
                                        ]})
