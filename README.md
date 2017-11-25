superpose3d
===========

##  Usage

```python
def Superpose3D(X_i     # <-- Nx3 array of coords for the "frozen" point cloud
                x_i,    # <-- Nx3 array of coords for the "mobile" point cloud
                w_i=[], # <-- optional weights for the calculation of RMSD
                        #     (default w_i = 1 for all i)
                allow_rescale=False)  #<--attempt to rescale mobile point cloud?
```

Superpose3D() takes two ordered lists (or numpy arrays) of xyz coordinates
(*of the same length*) representing points in a point cloud ("X_i", and "x_i"),
Treating them as rigid objects, "Superpose3D()" attempts to superimpose
them using **rotations**, **translations**, and (optionally) **scale**
transformations in order to minimize the root-mean-squared-distance (RMSD)
between corresponding points from either point cloud, where RMSD is defined as:

```
   RMSD =  sqrt( (Sum_i  w_i * |X_i - c*(R*x_i + T)|^2)  /  (Sum_j w_j) )
```
where:
```
   T = a translation vector (a numpy array containing x,y,z offsets),
   R = a rotation matrix    (a 3x3 numpy array),
   c = a scalar             (a number)
```
This function returns a 4-tuple containing the optimal values of:
```
   (RMSD, T, R, c)
```
This function implements a more general variant of the method from this paper:
R. Diamond, (1988)
"A Note on the Rotational Superposition Problem",
 Acta Cryst. A44, pp. 211-216
The version in the paper only considers rotation and translation and does not
allow the coordinates of either object to be scaled (multiplied by a scalar).


## Requirements

superpose3d depends on numpy

## License

superpose3d is available under the terms of the open-source 3-clause
BSD license.  (See `LICENSE.md`.)