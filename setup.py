from distutils.core import setup

setup(
    name='kmertools',
    version='0.1.0',
    author='Simone Longo',
    author_email='s.longo@utah.edu',
    packages=['kmertools'],
    license='LICENSE.txt',
    description='Tools for processing VCF files in parallel in addition to k-mer search and analysis',
    long_description=open('README.rst').read(),
    install_requires=[
        "pyfaidx",
        "cyvcf2",
        "pandas",
        "numpy",
    ],
)
