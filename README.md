[![Build Status](https://travis-ci.org/jewettaij/superpose3d.svg?branch=master)](./.travis.yml)
[![GitHub](https://img.shields.io/github/license/jewettaij/superpose3d)](./LICENSE.md)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/superpose3d)](https://pypistats.org/packages/superpose3d)
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/jewettaij/superpose3d)]()



superpose3d
===========

Note: There is a ***C++ version*** of this repository
[***here***](https://github.com/jewettaij/superpose3d_cpp).

##  Usage

```python
def Superpose3D(X,    # <-- Nx3 array of coords for the "frozen" point cloud
                x,    # <-- Nx3 array of coords for the "mobile" point cloud
                w = None, # <-- an optional array of N weights
                          #     (If w=None, equal weights will be used)
                allow_rescale=False,  #<--attempt to rescale mobile point cloud?
                q = None) # <-- optional: store the quaternion for rotation here
```

Superpose3D() takes two ordered lists (or numpy arrays) of xyz coordinates
(*of the same length*, **N**) representing points in a point cloud (**X** and
**x**). Treating them as rigid objects, "Superpose3D()" attempts to superimpose
them using **rotations**, **translations**, and (optionally) **scale**
transformations in order to minimize the root-mean-squared-distance (RMSD)
between corresponding points from either point cloud, where RMSD is defined as:

<img src="http://latex.codecogs.com/gif.latex?\large&space;RMSD=\sqrt{\,\frac{1}{N}\,\sum_{n=1}^N\,\,\sum_{i=1}^3 \left|X_{ni}-\left(\sum_{j=1}^3 cR_{ij}x_{nj}+T_i\right)\right|^2}"/>

...where:
```
   T_j  = a translation vector (a 1-D numpy array containing x,y,z displacements),
   R_ij = a rotation matrix    (a 3x3 numpy array whose determinant = 1),
    c   = a scalar             (a number, 1 by default)
```
This function returns a 4-tuple containing the optimal values of:
```
   (RMSD, T, R, c)
```

### Rotation angles, axes, and quaternions
If the rotation angle and axis are also needed, you can specify an
optional *q* argument (an array of size 4).  After invoking Superpose(),
*q* will store the
[quaternion corresponding to rotation *R*](https://en.wikipedia.org/wiki/Quatern
ions_and_spatial_rotation).  The first element of *q* will store *cos(θ/2)*
(where *θ* is the rotation angle).  The remaining 3 elements of *q* will store
the axis of rotation (with length *sin(θ/2)*).

### Weighted RMSD
A weighted version of the RMSD minimization algorithm is also available
if the caller supplies an extra argument specifying the weight of every
point in the cloud (*w<sub>n</sub>*).  In that case, RMSD is defined as:
<img src="http://latex.codecogs.com/gif.latex?\large&space;RMSD=\sqrt\left\sum_{n=1}^N\,w_n\,\sum_{i=1}^3 \left|X_{ni}-\left(\sum_{j=1}^3 c R_{ij}x_{nj}+T_i\right)\right|^2\quad\middle/\quad\sum_{n=1}^N w_n}\right}"/>

### Scale transformations
This function implements a more general variant of the method from this paper:
R. Diamond, (1988)
"A Note on the Rotational Superposition Problem",
 Acta Cryst. A44, pp. 211-216.

This version has been augmented slightly to support scale transformations.  (I.E. multiplication by scalars.  This can be useful for the registration of two different annotated volumetric 3-D images of the same object taken at different magnifications.)

Note that if you enable scale transformations (i.e. if *allow_rescale=True*), you should be wary if the function returns a negative **c** value.  Negative **c** values correspond to inversions (reflections).  For this reason, if you are using this function to compare the conformations of molecules, you should probably set *allow_rescale=False*.  This will prevent matching a molecule with its stereoisomer.

## Installation using pip

    pip install .
    pip install -r requirements.txt

Later, you can uninstall superpose3d using:

    pip uninstall superpose3d

## Requirements

superpose3d depends on numpy

## License

superpose3d is available under the terms of the [MIT license](LICENSE.md).
