# -*- coding: utf-8 -*-
#
import glob
import os
import codecs

from setuptools import setup, find_packages

# https://packaging.python.org/single_source_version/
base_dir = os.path.abspath(os.path.dirname(__file__))
about = {}
with open(os.path.join(base_dir, "quadpy", "__about__.py"), "rb") as f:
    exec(f.read(), about)


def read(fname):
    return codecs.open(os.path.join(base_dir, fname), encoding="utf-8").read()


directories = [
    os.path.join("quadpy", "triangle", "papanicolopulos"),
    os.path.join("quadpy", "triangle", "vioreanu_rokhlin"),
    os.path.join("quadpy", "triangle", "witherden_vincent"),
    os.path.join("quadpy", "triangle", "xiao_gimbutas"),
    #
    os.path.join("quadpy", "sphere", "lebedev"),
    #
    os.path.join("quadpy", "tetrahedron", "xiao_gimbutas"),
]

data_files = [
    (directory, glob.glob(os.path.join(directory, "*.json")))
    for directory in directories
]

setup(
    name="quadpy",
    version=about["__version__"],
    packages=find_packages(),
    url="https://github.com/nschloe/quadpy",
    author=about["__author__"],
    author_email=about["__email__"],
    install_requires=[
        "matplotlib",
        "numpy",
        "orthopy >=0.5, <0.6",
        "pipdate >=0.2.0, <0.3.0",
        "scipy",
        "sympy",
    ],
    description="Numerical integration, quadrature for various domains",
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    license=about["__license__"],
    classifiers=[
        about["__license__"],
        about["__status__"],
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    data_files=data_files,
)
