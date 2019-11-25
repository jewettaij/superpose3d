#!/usr/bin/env python

import numpy as np
from math import *
from superpose3d import Superpose3D

def test_superpose3d():
    X=[[0,0,0],[0,0,1],[0,1,1],[1,1,1]]
    x=[[0,0,0],[0,1,0],[0,1,-1],[1,1,-1]]  # (=X after rotation around the Z axis)

    Xshifted = [ [X[i][0],X[i][1]+100, X[i][2]] for i in range(0,len(X))]
    xscaled  = [ [2*x[i][0],    2*x[i][1],2*x[i][2]] for i in range(0,len(x))]
    xscshift = [ [2*x[i][0]+200,2*x[i][1],2*x[i][2]] for i in range(0,len(x))]

    result = Superpose3D(X,x)
    assert(abs(result[0]) < 1.0e-6)



    result = Superpose3D(X,xscshift, None, True)

    # Does the RMSD returned in result[0] match the RMSD calculated manually?
    R = np.matrix(result[1])              # rotation matrix
    T = np.matrix(result[2]).transpose()  # translation vector (3x1 matrix)
    c = result[3]                         # scalar
    _x = np.matrix(xscshift).transpose()
    _xprime = c*R*_x + T
    xprime = np.array(_xprime.transpose()) # convert to length 3 numpy array
    RMSD = 0.0

    #print(X)
    #print(xprime)

    for i in range(0, len(X)):
        RMSD += ((X[i][0] - xprime[i][0])**2 +
                 (X[i][1] - xprime[i][1])**2 +
                 (X[i][2] - xprime[i][2])**2)

    RMSD = sqrt(RMSD / len(X))
    assert(abs(RMSD - result[0]) < 1.0e-6)


def main():
    test_superpose3d()


if __name__ == '__main__':
    main()
