from setuptools import setup, find_packages

import ta2_hrr_2019 # This should trigger an error if paths.py is missing

setup(name='ta2_hrr_2019',
      version='0.0.0',
      packages=find_packages(),
      zip_safe=False,
      include_package_data=True,
      )
