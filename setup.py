import sys
import glob
from setuptools import setup

setup(name='PMF3D',
      version = '1.0',
      description = 'Protein-ligand binding free energy maps', 
      packages=['PMF3D',
                'PMF3D.scripts'],
      package_dir={'PMF3D':'lib',
                   'PMF3D.scripts':'scripts'},
      install_requires=['ipython==0.13'],
      scripts=filter(lambda elem: '_' not in elem, glob.glob('scripts/*')))
