from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
    name='highlighting_atoms',
    version='0.5',
    description='Highlighting atoms function',
    license='GNU',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Pablo Rodríguez Belenguer',
    author_email='parodbe@gmail.com',
    url='https://github.com/phi-grib/highlighting_atoms',
    download_url='https://github.com/phi-grib/highlighting_atoms',
    packages=find_packages())
