from setuptools import setup

setup(
  
  name='superpose3d',

    packages=['superpose3d'],

  description='Diamond\'s 1988 rotational superposition algorithm (+scale tranforms)',

  author='Andrew Jewett',

  author_email='jewett.aij@gmail.com',

  url='https://github.com/jewettaij/superpose3d',

  download_url='https://github.com/jewettaij/superpose3d/archive/v0.0.3.zip',

  version='0.0.3',

  # BSD 3-Clause License:
  # - http://choosealicense.com/licenses/bsd-3-clause
  # - http://opensource.org/licenses/BSD-3-Clause

  license='BSD',

  classifiers=['Development Status :: 3 - Alpha',
               'License :: OSI Approved :: BSD License',
               'Environment :: Console',
               'Operating System :: MacOS :: MacOS X',
               'Operating System :: POSIX :: Linux',
               'Operating System :: Microsoft :: Windows',
               'Programming Language :: Python',
               'Topic :: Scientific/Engineering'],

  zip_safe=True,
  include_package_data=True
)
