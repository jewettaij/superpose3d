[![Build Status](https://travis-ci.org/jewettaij/superpose3d.svg?branch=master)](./.travis.yml)
[![GitHub](https://img.shields.io/github/license/jewettaij/superpose3d)](./LICENSE.md)
[![PyPI - Downloads](https://img.shields.io/pypi/dm/superpose3d)](https://pypistats.org/packages/moltemplate)
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/jewettaij/superpose3d)]()



superpose3d
===========

*Note: C++ version [here](https://github.com/jewettaij/superpose3d_cpp)*

##  Usage

```python
def Superpose3D(X_i,    # <-- Nx3 array of coords for the "frozen" point cloud
                x_i,    # <-- Nx3 array of coords for the "mobile" point cloud
                w_i=None, #<- optional weights for the calculation of RMSD
                          #   (default w_i = 1 for all i)
                allow_rescale=False)  #<--attempt to rescale mobile point cloud?
```

Superpose3D() takes two ordered lists (or numpy arrays) of xyz coordinates
(*of the same length*, **N**) representing points in a point cloud (**X_i** and
**x_i**). Treating them as rigid objects, "Superpose3D()" attempts to superimpose
them using **rotations**, **translations**, and (optionally) **scale**
transformations in order to minimize the root-mean-squared-distance (RMSD)
between corresponding points from either point cloud, where RMSD is defined as:

<img src="http://latex.codecogs.com/gif.latex?\large&space;RMSD=\left(\frac{\sum_{i=1}^n\,w_i\,|X_i-\sum_{j=1}^n(cR_{ij}x_j+T_i)|^2}{\sum_{i=1}^nw_i}\right)^{\frac{1}{2}}"/>

If *w<sub>i</sub>* are omitted (ie. if *w<sub>i</sub> = nullptr*),
then equal weights are used.  In that case:

<img src="http://latex.codecogs.com/gif.latex?\large&space;RMSD=\left(\frac{1}{n}\,\sum_{i=1}^n\,|X_i-\sum_{j=1}^n (cR_{ij}x_j+T_i)|^2\right)^{\frac{1}{2}}"/>

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

## Installation using pip

    pip install .
    pip install -r requirements.txt

Later, you can uninstall superpose3d using:

    pip uninstall superpose3d

## Requirements

superpose3d depends on numpy

## License

superpose3d is available under the terms of the [MIT license](LICENSE.md).
