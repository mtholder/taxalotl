#!/usr/bin/env python
import sys
import os

# setup.py largely based on
#   http://hynek.me/articles/sharing-your-labor-of-love-pypi-quick-and-dirty/
# Also see Jeet Sukumaran's DendroPy

###############################################################################
# setuptools/distutils/etc. import and configuration
try:
    # noinspection PyUnresolvedReferences
    import ez_setup

    try:
        ez_setup_path = " ('" + os.path.abspath(ez_setup.__file__) + "')"
    except OSError:
        ez_setup_path = ""
    sys.stderr.write("using ez_setup%s\n" % ez_setup_path)
    ez_setup.use_setuptools()
    # noinspection PyUnresolvedReferences
    import setuptools

    try:
        setuptools_path = " ('" + os.path.abspath(setuptools.__file__) + "')"
    except OSError:
        setuptools_path = ""
    sys.stderr.write("using setuptools%s\n" % setuptools_path)
    # noinspection PyUnresolvedReferences
    from setuptools import setup, find_packages
except ImportError as e:
    sys.stderr.write("using distutils\n")
    from distutils.core import setup

    sys.stderr.write("using canned package list\n")
    PACKAGES = ['taxalotl',
                ]
    EXTRA_KWARGS = {}
else:
    sys.stderr.write("searching for packages\n")
    PACKAGES = find_packages()
    # noinspection PyTypeChecker
    EXTRA_KWARGS = dict(
        include_package_data=True,
        test_suite="taxalotl.test"
    )

EXTRA_KWARGS["zip_safe"] = False
ENTRY_POINTS = {}
EXTRA_KWARGS['scripts'] = ['taxalotl-cli.py',
                           ]
setup(
    name='taxalotl',
    version='0.0dev',  # sync with __version__ in peyotl/__init__.py
    description='Taxonomy merging tool',
    long_description=(open('README.md').read()),
    url='https://github.com/mtholder/taxalotl',
    license='BSD',
    author='Mark T. Holder',
    py_modules=['taxalotl'],
    install_requires=['setuptools',
                      'requests>=2.2.1', 'ez_setup',
                      ],
    packages=PACKAGES,
    entry_points=ENTRY_POINTS,
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    **EXTRA_KWARGS
)
