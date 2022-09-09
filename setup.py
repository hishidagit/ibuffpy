from setuptools import find_packages, setup

setup(
    name='SSA',
    packages=find_packages('src'),
    package_dir={"": "src"},
)