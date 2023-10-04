from setuptools import find_packages, setup
# import ssapy


INSTALL_REQUIRES = [
    'numpy>=1.23.5',
    'scipy>=1.10.1',
    'matplotlib>=3.7.3',
    'pandas>=1.5.3',
    'tqdm>=4.66.1',
    'sympy>=1.12',
    'networkx>=3.1',
    'cobra!=0.26.0',
    'scikit-learn>=1.3.1',
]

with open('./README.md', 'r') as f:
    long_description = f.read()

setup(
    # version=ssapy.__version__,
    version='0.0.3',
    name='ssapy',
    description='A Python package for appling Structural Sensitivity Analysis to reaction networks.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url = 'https://github.com/hishidagit/SSApy',
    author='Atsuki Hishida',
    author_email='ahishida18@gmail.com',
    license='MIT',
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"],
    project_urls={
        "Homepage" : "https://github.com/hishidagit/SSApy",
        "Bug Tracker" : "https://github.com/hishidagit/SSApy/issues"},
    packages=find_packages('src'),
    package_dir={"": "src"},
    install_requires=INSTALL_REQUIRES,
    python_requires='>=3.8',
)