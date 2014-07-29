import sys
import glob
from setuptools import setup

setup(name='PMF3D',
      author='Morgan Lawrenz', 
      author_email='morganlawrenz@gmail.com',
      version = '1.0',
      description = 'Protein-ligand binding 3-D free energy maps', 
      py_modules = ['PMF3D'],
      scripts=glob.glob('RunPMF.py'), )
