# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='pyHeatBJ',
    version='0.1.0',
    description='Sample package for Python',
    long_description=readme,
    author='Julien Horsin & Berenice Sinopoli--Pal',
    author_email='julien.horsin@mines-paristech.fr',
    url='https://github.com/BereniceSinopoliPal/Calcul_Molonari_Julien_Berenice',
    license=license,
    packages=find_packages()
)

