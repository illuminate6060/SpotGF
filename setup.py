from setuptools import setup, find_packages
import os

requires = ['numpy',
            'pot==0.8.2',
            'pandas',
            'matplotlib',
            'scipy',
            'descartes',
            'alphashape',
            'shapely'
            ]

setup(
      name='GeneFilter', 
      version='1.0', 
      author='Lin Du',
      author_email='3051065449@qq.com',
      license="MIT",
      description="GF for de-noising spatially resolved transcriptomics data based on optimal transport theory",
      # long_description=open('README').read(),
      packages=find_packages(), 
      install_requires = requires,
      python_requires = '>=3'
)