#!/usr/bin/env python

from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()

setup(
    name="cine",
    version="0.2",
    author="Miguel de Val-Borro",
    author_email="miguel.deval@gmail.com",
    url="https://github.com/migueldvb/cine",
    packages=["cine"],
    install_requires=[
        'numpy',
    ],
    scripts=['bin/cine'],
    description="Calculate infrared pumping rates by solar radiation",
    long_description=readme(),
    package_data={"": ["LICENSE"]},
    include_package_data=True,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    test_suite='nose.collector',
    tests_require=['nose'],
)
