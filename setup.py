import codecs
import os

from setuptools import find_packages, setup

# https://packaging.python.org/single_source_version/
base_dir = os.path.abspath(os.path.dirname(__file__))
about = {}
with open(os.path.join(base_dir, "quadpy", "__about__.py"), "rb") as f:
    exec(f.read(), about)


def read(fname):
    return codecs.open(os.path.join(base_dir, fname), encoding="utf-8").read()


setup(
    name="quadpy",
    version=about["__version__"],
    packages=find_packages(),
    url="https://github.com/nschloe/quadpy",
    author=about["__author__"],
    author_email=about["__email__"],
    install_requires=["numpy", "orthopy >=0.5, <0.6", "scipy", "sympy"],
    extras_require={"all": ["matplotlib"], "plot": ["matplotlib"]},
    description="Numerical integration, quadrature for various domains",
    long_description=read("README.md"),
    long_description_content_type="text/markdown",
    license=about["__license__"],
    classifiers=[
        about["__license__"],
        about["__status__"],
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    package_data={"": ["*.json"]},
)
