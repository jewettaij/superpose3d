# -*- coding: utf-8 -*-

from math import *
import numpy as np
from numpy import linalg as LA

def Superpose3D(aaXf_orig,   # <-- coordinates for the "frozen" object
                aaXm_orig,   # <-- coordinates for the "mobile" object
                aWeights=None, #<- optional weights for the calculation of RMSD
                allow_rescale=False): # <-- attempt to rescale mobile object?
    """
    Superpose3D() takes two lists of xyz coordinates, (of the same length)
    and attempts to superimpose them using rotations, translations, and 
    (optionally) rescale operations in order to minimize the 
    root-mean-squared-distance (RMSD) between them.  
    These operations should be applied to the "aaXm_orig" argument.
    This function returns a tuple containing:
      (RMSD, optimal_translation, optimal_rotation, and optimal_scale_factor)
    This function implements a more general variant of the method from:
    R. Diamond, (1988)
    "A Note on the Rotational Superposition Problem", 
    Acta Cryst. A44, pp. 211-216
    This version has been augmented slightly.  The version in the original 
    paper only considers rotation and translation and does not allow the 
    coordinates of either object to be rescaled (multiplication by a scalar).

    (Additional documentation can be found at
     https://pypi.org/project/superpose3d/ )

    """

    assert(len(aaXf_orig) == len(aaXm_orig))
    N = len(aaXf_orig)
    if (aWeights == None) or (len(aWeights) == 0):
        aWeights = np.full(N,1.0)

    # Find the center of mass of each object:
    aCenter_f = np.zeros(3)
    aCenter_m = np.zeros(3)
    sum_weights = 0.0
    for n in range(0, N):
        for d in range(0, 3):
            aCenter_f[d] += aaXf_orig[n][d]*aWeights[n]
            aCenter_m[d] += aaXm_orig[n][d]*aWeights[n]
        sum_weights += aWeights[n]
    for d in range(0, 3):
        aCenter_f[d] /= sum_weights
        aCenter_m[d] /= sum_weights

    # Subtract the centers-of-mass from the original coordinates for each object
    aaXf = np.empty((N,3))
    aaXm = np.empty((N,3))
    aaXf[0][0] = 0.0
    for n in range(0, N):
        for d in range(0, 3):
            aaXf[n][d] = aaXf_orig[n][d] - aCenter_f[d]
            aaXm[n][d] = aaXm_orig[n][d] - aCenter_m[d]

    Rgf = 0.0  # <-- the RMS size of the particles in the frozen object aaXf
    Rgm = 0.0  # <-- the RMS size of the particles in the mobile object aaXm

    if allow_rescale:
        # Optional: For numerical stability, we might as well rescale the
        # coordinates initially to make sure they have the same approximate
        # scale before we attempt to superimpose them.
        # This is only necessary if one object is much bigger than the other
        # (ie. by several orders of magnitude).
        # Note: This is NOT the optimal scale factor.
        #       (That must be determined later.)
        for n in range(0, N):
            for d in range(0, 3):
                Rgf += aWeights[n]*((aaXf[n][d])**2)
                Rgm += aWeights[n]*((aaXm[n][d])**2)
        Rgf = sqrt(Rgf / sum_weights)
        Rgm = sqrt(Rgm / sum_weights)

        for n in range(0, N):
            for d in range(0, 3):
                aaXf[n][d] /= Rgf
                aaXm[n][d] /= Rgm

    # Calculate the "M" array from the Diamond paper (equation 16)
    M = np.zeros((3,3))
    for n in range(0, N):
        for i in range(0, 3):
            for j in range(0, 3):
                M[i][j] += aWeights[n] * aaXm[n][i] * aaXf[n][j]

    # Calculate Q (equation 17)
    traceM = 0.0
    for i in range(0, 3):
        traceM += M[i][i]

    Q = np.empty((3,3))
    for i in range(0, 3):
        for j in range(0, 3):
            Q[i][j] = M[i][j] + M[j][i]
            if i==j:
                Q[i][j] -= 2.0 * traceM

    # Calculate V (equation 18)
    V = np.empty(3)
    V[0] = M[1][2] - M[2][1];
    V[1] = M[2][0] - M[0][2];
    V[2] = M[0][1] - M[1][0];

    # Calculate "P" (equation 22)
    P = np.empty((4,4))
    for i in range(0,3):
        for j in range(0,3):
            P[i][j] = Q[i][j]
    P[0][3] = V[0]
    P[3][0] = V[0]
    P[1][3] = V[1]
    P[3][1] = V[1]
    P[2][3] = V[2]
    P[3][2] = V[2]
    P[3][3] = 0.0

    aEigenvals, aaEigenvects = LA.eigh(P)

    #http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.eigh.html

    eval_max = aEigenvals[0]
    i_eval_max = 0
    for i in range(1,4):
        if aEigenvals[i] > eval_max:
            eval_max = aEigenvals[i]
            i_eval_max = i

    # The vector "p" contains the optimal rotation (in quaternion format)
    p = np.empty(4)
    p[0] = aaEigenvects[0][ i_eval_max ]
    p[1] = aaEigenvects[1][ i_eval_max ]
    p[2] = aaEigenvects[2][ i_eval_max ]
    p[3] = aaEigenvects[3][ i_eval_max ]

    # normalize the vector
    # (It should be normalized already, but just in case it is not, do it again)
    pnorm = np.linalg.norm(p)
    for i in range(0,4):
        p[i] /= pnorm

    # Finally, calculate the rotation matrix corresponding to "p"
    # (convert a quaternion into a 3x3 rotation matrix)
    aaRotate = np.empty((3,3))

    aaRotate[0][0] =  (p[0]*p[0])-(p[1]*p[1])-(p[2]*p[2])+(p[3]*p[3])
    aaRotate[1][1] = -(p[0]*p[0])+(p[1]*p[1])-(p[2]*p[2])+(p[3]*p[3])
    aaRotate[2][2] = -(p[0]*p[0])-(p[1]*p[1])+(p[2]*p[2])+(p[3]*p[3])
    aaRotate[0][1] = 2*(p[0]*p[1] - p[2]*p[3]);
    aaRotate[1][0] = 2*(p[0]*p[1] + p[2]*p[3]);
    aaRotate[1][2] = 2*(p[1]*p[2] - p[0]*p[3]);
    aaRotate[2][1] = 2*(p[1]*p[2] + p[0]*p[3]);
    aaRotate[0][2] = 2*(p[0]*p[2] + p[1]*p[3]);
    aaRotate[2][0] = 2*(p[0]*p[2] - p[1]*p[3]);
    
    pPp = eval_max

    # Optional: Decide the scale factor, c
    c = 1.0   # by default, don't rescale the coordinates
    if allow_rescale:
        Waxaixai = 0.0
        WaxaiXai = 0.0
        for a in range(0, N):
            for i in range(0, 3):
                Waxaixai += aWeights[a] * aaXm[a][i] * aaXm[a][i]
                WaxaiXai += aWeights[a] * aaXm[a][i] * aaXf[a][i]
        c = (WaxaiXai + pPp) / Waxaixai
        # Recall that we previously divided the two sets of coordinates by Rgm
        # and Rgf respectively. (I thought it might improve numerical stability)
        # Before returning "c" to the caller, we need to incorporate those
        # factors into "c" as well.
        c *= Rgf / Rgm
        pPp *= Rgf * Rgm
        # And, lastly, undo this before calculating E0 below
        for n in range(0, N):
            for d in range(0, 3):
                aaXf[n][d] *= Rgf
                aaXm[n][d] *= Rgm
    # Finally compute the RMSD between the two coordinate sets:
    # First compute E0 from equation 24 of the paper
    E0 = 0.0
    for n in range(0, N):
        for d in range(0, 3):
            # (remember to include the scale factor "c" that we inserted)
            E0 += aWeights[n] * ((aaXf[n][d] - c*aaXm[n][d])**2)
    sum_sqr_dist = E0 - c*2.0*pPp
    if sum_sqr_dist < 0.0:
        sum_sqr_dist = 0.0
    rmsd = sqrt(sum_sqr_dist/sum_weights)

    # Lastly, calculate the translational offset:
    # Recall that:
    #RMSD=sqrt((Σ_i  w_i * |X_i - (Σ_j c*R_ij*x_j + T_i))|^2) / (Σ_j w_j))
    #    =sqrt((Σ_i  w_i * |X_i - x_i'|^2) / (Σ_j w_j))
    #  where
    # x_i' = Σ_j c*R_ij*x_j + T_i
    #      = Xcm_i + c*R_ij*(x_j - xcm_j)
    #  and Xcm and xcm = center_of_mass for the frozen and mobile point clouds
    #                  = aCenter_f[]       and       aCenter_m[],  respectively
    # Hence:
    #  T_i = Xcm_i - Σ_j c*R_ij*xcm_j  =  aTranslate[i]

    aTranslate = np.empty(3)
    for i in range(0,3):
        aTranslate[i] = aCenter_f[i]
        for j in range(0,3):
            aTranslate[i] -= c*aaRotate[i][j]*aCenter_m[j]

    # An alternate method to compute "aTranslate" using numpy matrices:
    #Rmatrix = np.matrix(aaRotate)
    #TcolumnVec = np.matrix(np.empty((3,1))) # 3x1 numpy matrix<->[[0],[0],[0]]
    #for d in range(0,3):
    #    TcolumnVec[d][0] = -aCenter_m[d]
    #TcolumnVec = c * Rmatrix * TcolumnVec
    #for d in range(0,3):
    #    TcolumnVec[d][0] += aCenter_f[d]
    # #Turn the column vector back into an ordinary numpy array of size 3:
    #aTranslate = np.array(TcolumnVec.transpose())[0]

    return rmsd, aaRotate, aTranslate, c

