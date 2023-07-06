from setuptools import setup, find_packages
import os

requires = ['numpy',
            'pot==0.8.2',
            'pandas',
            'matplotlib',
            'scipy',
            'descartes',
            ]

setup(
      name='GeneFilter', 
      version='0.1', 
      author='dulin',
      author_email='dulin@genomics.cn',
      license="MIT",
      description="Optimal Transport Method-Based Gene Filter (GF) Denoising Algorithm for Enhancing Spatially Resolved Transcriptomics Data",
      # long_description=open('README').read(),
      packages=find_packages(), 
      install_requires = requires,
      python_requires = '>=3'
)