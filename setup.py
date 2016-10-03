# -*- coding: utf-8 -*-
#
from distutils.core import setup
import os
import codecs

from quadrature import __version__, __license__, __author__, __email__


def read(fname):
    try:
        content = codecs.open(
            os.path.join(os.path.dirname(__file__), fname),
            encoding='utf-8'
            ).read()
    except Exception:
        content = ''
    return content

setup(
    name='quadrature',
    version=__version__,
    packages=['quadrature'],
    url='https://github.com/nschloe/quadrature',
    download_url='https://pypi.python.org/pypi/quadrature',
    author=__author__,
    author_email=__email__,
    requires=['numpy'],
    description='numerical integration schemes',
    long_description=read('README.rst'),
    license=__license__,
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics'
        ]
    )
