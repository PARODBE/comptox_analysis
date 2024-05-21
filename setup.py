from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
    name='Highlighting_atoms_rdkit',
    version='0.1',
    description='Highlighting atoms function',
    license='GNU',
#    long_description=long_description,
    author='Pablo Rodríguez Belenguer',
    author_email='parodbe@gmail.com',
    url='https://github.com/phi-grib/Highlighting_atoms_rdkit',
    download_url='https://github.com/phi-grib/Highlighting_atoms_rdkit.git',
    packages=find_packages(),
    # If any package contains *.txt or *.rst files, include them:
    # package_data={'': ['*.yaml', '*.yml']},
)
