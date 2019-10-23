from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name = 'bucoffea',
    version = '0.0.1',
    url = 'https://github.com/bu-cms/bucoffea',
    author = 'Andreas Albert',
    author_email = 'andreas.albert@cern.ch',
    description = 'Analysis using Coffea on NanoAOD',
    packages = find_packages(),    
    install_requires = requirements,
    scripts=[
        'bucoffea/execute/buexec',
        'bucoffea/execute/bumon',
        'bucoffea/scripts/bumerge'
        ],
)
