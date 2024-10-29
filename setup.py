from setuptools import setup, find_packages
import os

requires = ['numpy',
            'pot==0.8.2',
            'pandas',
            'matplotlib',
            'scipy',
            'descartes',
            'alphashape',
            'shapely<2.0.0'
            ]

setup(
      name='SpotGF', 
      version='1.0', 
      author='Lin Du',
      author_email='3051065449@qq.com',
      license="MIT",
      description="SpotGF: Denoising Spatially Resolved Transcriptomics Data Using an Optimal Transport-Based Gene Filtering Algorithm",
      long_description=open('README.md').read(),
      packages=find_packages(), 
      install_requires = requires,
      python_requires = '>=3'
)
