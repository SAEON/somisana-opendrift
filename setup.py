from setuptools import setup, find_packages
# bare minimum for being able to install our functions in a conda environment - needed when doing
# pip install --no-deps -e . 
setup(
    name='somisana-opendrift',
    version='0.1',
    packages=find_packages(),
)

