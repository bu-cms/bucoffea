from setuptools import setup, find_packages

setup(
    name = 'bucoffea',
    version = '0.0.1',
    url = 'https://github.com/AndreasAlbert/monojetcoffea.git',
    author = 'Andreas Albert',
    author_email = 'andreas.albert@cern.ch',
    description = 'Monojet analysis using Coffea on NanoAOD',
    packages = find_packages(),    
    install_requires = ['coffea'],
)