#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import setuptools

from setuptools import setup, Extension

setup(name = 'nlg_pow',
    version = '1.0.0',
    description = 'Computation of analogical power',
    author = 'Yves Lepage',
    author_email = 'yves.lepage@waseda.jp',
    ext_modules =	[Extension('_nlg_pow',
                        sources = ['nlg_pow.i', 'nlg_pow.c'],
                        extra_compile_args = ['-std=c99'],
                    )],
    packages = setuptools.find_packages(),
    python_requires = '>=3.11'
    )
