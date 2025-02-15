from setuptools import setup, find_packages
from os import path
working_directory = path.abspath(path.dirname(__file__))

with open(path.join(working_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='c4_open',
    version='0.1.3',
    author='Carsten Hennig',
    author_email='carsten_hennig@gmx.de',
    description='A package for current carrying capacity calculation of cables',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/carstenhennigberlin/c4_open",
    packages=find_packages(where='c4_open'),
    package_dir={'': 'c4_open'},
    include_package_data=True,
    package_data={
        'c4_open': ['data/*.db'],
    },
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
        'ezdxf',
        'openpyxl',
        'xlrd',
        'scipy',
    ]
)