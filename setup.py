#!/usr/bin/env python3

from setuptools import setup

setup(
  name="helper_functions",
  version='0.1.0',
  description="Helper functions for genomics analysis",
  long_description=open('README.md').read(),
  long_description_content_type='text/markdown',
  url='https://github.com/bkinnersley/misc',
  author="Ben Kinnersley",
  license="GPL",
  python_requires='>= 3.7',
  author_email="b.kinnersley@ucl.ac.uk",
  packages=["functions"],
  install_requires=["datetime"],
  include_package_data=True,
  zip_safe=False
)
