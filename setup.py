#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

__author__ = 'Yuki Kuriya'

setup(
    name='THP_pro_pathway_select',  
    version='1.0',  
    description='An example of system analysis and pathway selection by Monte calro simulation',  
    author='Yuki Kuriya', 
    author_email='yukuriya3.kobe@gmail.com',  
    url='https://github.com/yukuriya3/THP_pro_pathway_select', 
    classifiers=[ 
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3 :: Only',
    ],
    packages=find_packages(exclude=('tests', 'docs')), 
    include_package_data=True,  
    keywords=['Monte Carlo simulation'], 
    license='MIT License', 
    install_requires=[ # package dependency 
        'numpy',
        'scipy',
        'matplotlib',
        'pandas',
        'csv',
    ],
)

