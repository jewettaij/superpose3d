# -*- coding: utf-8 -*-

from setuptools import setup

setup(
  
  name='superpose3d',

  packages=['superpose3d'],

  description='Diamond\'s 1988 rotational superposition algorithm (+scale tranforms)',
  long_description='''Register 3-D point clouds using rotation, translation, and scale transformations.

##  Usage

```python
def Superpose3D(X,    # <-- Nx3 array of coords for the "frozen" point cloud
                x,    # <-- Nx3 array of coords for the "mobile" point cloud
                w = None, #<- optional weights for the calculation of RMSD
                          #   (default w[n] = 1 for all n)
                allow_rescale=False)  #<--attempt to rescale mobile point cloud?
```

Superpose3D() takes two ordered lists (or numpy arrays) of xyz coordinates
(*of the same length*, **N**) representing points in a point cloud (**X** and
**x**). Treating them as rigid objects, "Superpose3D()" attempts to superimpose
them using **rotations**, **translations**, and (optionally) **scale**
transformations in order to minimize the root-mean-squared-distance (RMSD)
between corresponding points from either point cloud, where RMSD is defined as:

<img src="http://latex.codecogs.com/gif.latex?\large&space;RMSD=\sqrt\left\sum_{n=1}^N\,w_n\,\sum_{i=1}^3 \left|X_{ni}-\left(\sum_{j=1}^3 c R_{ij}x_{nj}+T_i\right)\right|^2\quad\middle/\quad\sum_{n=1}^N w_n}\right}"/>

If *w<sub>n</sub>* are omitted (ie. if *w<sub>n</sub> = None*),
then equal weights are used.  In that case:

<img src="http://latex.codecogs.com/gif.latex?\large&space;RMSD=\sqrt{\,\frac{1}{n}\,\sum_{n=1}^N\,\,\sum_{i=1}^3 \left|X_{ni}-\left(\sum_{j=1}^3 cR_{ij}x_{nj}+T_i\right)\right|^2}"/>

...where:

```
   T_j  = a translation vector (a 1-D numpy array containing x,y,z displacements),
   R_ij = a rotation matrix    (a 3x3 numpy array whose determinant = 1),
    c   = a scalar             (a number)
```
This function returns a 4-tuple containing the optimal values of:
```
   (RMSD, T, R, c)
```
This function implements a more general variant of the method from this paper:
R. Diamond, (1988)
"A Note on the Rotational Superposition Problem",
 Acta Cryst. A44, pp. 211-216.

This version has been augmented slightly to support scale transformations.  (I.E. multiplication by scalars.  This can be useful for the registration of two different annotated volumetric 3-D images of the same object taken at different magnifications.)

Note that if you enable scale transformations (i.e. if *allow_rescale=True*), you should be wary if the function returns a negative **c** value.  Negative **c** values correspond to inversions (reflections).  For this reason, if you are using this function to compare the conformations of molecules, you should probably set *allow_rescale=False*.  This will prevent matching a molecule with its stereoisomer.

Note: A C++ version of this repository is available at
https://github.com/jewettaij/superpose3d_cpp
''',

  long_description_content_type='text/markdown',

  author='Andrew Jewett',

  author_email='jewett.aij@gmail.com',

  url='https://github.com/jewettaij/superpose3d',

  download_url='https://github.com/jewettaij/superpose3d/archive/v0.2.10.zip',

  version='0.2.10',

  install_requires=[
      'numpy',
  ],
    
  keywords=['registration', '3d', 'structure-comparison', 'molecular-structure',
            'clem'],

  license='MIT',

  classifiers=['Development Status :: 4 - Beta',
               'License :: OSI Approved :: MIT License',
               'Environment :: Console',
               'Operating System :: MacOS :: MacOS X',
               'Operating System :: POSIX :: Linux',
               'Operating System :: Microsoft :: Windows',
               'Programming Language :: Python',
               'Topic :: Scientific/Engineering'],

  zip_safe=True,
  include_package_data=True
)
