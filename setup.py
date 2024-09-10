from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
    name='comptox_analysis',
    version='0.8',
    description='ML analysis',
    license='GNU',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Pablo Rodríguez Belenguer',
    author_email='parodbe@gmail.com',
    url='https://github.com/phi-grib/comptox_analysis',
    download_url='https://github.com/phi-grib/comptox_analysis',
    packages=find_packages())
