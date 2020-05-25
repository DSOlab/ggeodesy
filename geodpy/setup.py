from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))

setup(
  name='geodpy',
  version='1.0.1',  # Required
  description='A sample Python project',  # Optional
  url='https://github.com/pypa/sampleproject',  # Optional
  author='The Python Packaging Authority',  # Optional
  author_email='pypa-dev@googlegroups.com',  # Optional
  keywords='sample setuptools development',  # Optional
#  package_dir={'': 'geodpy'},  # Optional
#  packages=find_packages(where='geodpy'),  # Required
  packages=['geodpy', 'geodpy.ellipsoid', 'geodpy.units'],
  python_requires='>=3.5, <4',
  install_requires=[],  # Optional
)
# import setuptools
# 
# setuptools.setup(name='geodpy',
#     version='1.0',
#     description='Python Geodesy Library',
#     url='',
#     author='Xanthos Papanikolaou, Demitris Anastasiou',
#     author_email='xanthos@mail.ntua.gr, danast@mail.ntua.gr',
#     packages=['geodpy'],
#     install_requires=[]
#)
