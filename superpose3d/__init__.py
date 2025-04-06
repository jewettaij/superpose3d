# -*- coding: utf-8 -*-
"""
This module contains the definition of the Superpose3D() function used for 
registering two point clouds of known correspondence. (It is documented below.)

Note: The original version of this code contained for-loops.  Newer versions of
      this code use numpy expressions to avoid for-loops, however the original
      code remains in quoted comments because some users find it easier to read.
"""

import numpy as np
from numpy import linalg as LA
from numpy.typing import NDArray
from typing import Optional, Sequence, Tuple
#from math import *

def Superpose3D(aaXf_original: Sequence[float],   # <-- coordinates for the "frozen" object
                aaXm_original: Sequence[float],   # <-- coordinates for the "mobile" object
                # ---- optional arguments: ----
                aWeights_original: Optional[Sequence[float]] = None,  # optional weights for the calculation of RMSD
                allow_rescale: bool = False,      # attempt to rescale mobile point cloud?
                report_quaternion: bool = False,  # report rotation angle and axis?
    ) -> Tuple[float, NDArray, NDArray, float]:
    """
    Superpose3D() takes two lists of xyz coordinates, (of the same length)
    and attempts to superimpose them using rotations, translations, and 
    (optionally) rescale operations in order to minimize the 
    root-mean-squared-distance (RMSD) between them.  
    These operations should be applied to the "aaXm_orig" argument.
    This function returns a tuple containing:
      (RMSD, optimal_translation, optimal_rotation, and optimal_scale_factor)
    More detailed documentation can be found in the repository's README.md file.
    This function implements a more general variant of the method
    described in this paper R. Diamond, (1988) Acta Cryst. A44, pp. 211-216
      "A note on the rotational superposition problem"
      https://doi.org/10.1107/S0108767387010535
    This version has been augmented slightly.  The version in the original 
    paper only considers rotation and translation and does not allow the 
    coordinates of either object to be rescaled (multiplication by a scalar).
    """

    # Convert input lists as to numpy arrays

    aaXf_orig = np.array(aaXf_original)
    aaXm_orig = np.array(aaXm_original)

    if aaXf_orig.shape[0] != aaXm_orig.shape[0]:
        raise ValueError ("Inputs should have the same size.")

    N = aaXf_orig.shape[0]

    # Find the center of mass of each object.
    # ...new code (avoiding for-loops)
    # First convert weights into an array of the correct shape
    aWeights = np.full((N, 1), 1.0)
    if (aWeights_original is not None) and (len(aWeights_original) != 0):
        aWeights = np.array(aWeights_original).reshape(N, 1)
    aCenter_f: NDArray = np.sum(aaXf_orig * aWeights, axis=0)
    aCenter_m: NDArray = np.sum(aaXm_orig * aWeights, axis=0)
    sum_weights: float = np.sum(aWeights)
    """ # For reference, here's the same code using for-loops
    if (aWeights == None) or (len(aWeights) == 0):
        aWeights = np.full(N, 1.0)
    aCenter_f = np.zeros(3)
    aCenter_m = np.zeros(3)
    sum_weights = 0.0
    for n in range(0, N):
        for d in range(0, 3):
            aCenter_f[d] += aaXf_orig[n,d]*aWeights[n]
            aCenter_m[d] += aaXm_orig[n,d]*aWeights[n]
        sum_weights += aWeights[n]
    """

    # Rescale by 1/sum_weights (compensate multiplying by the weights earlier)
    if sum_weights != 0:
        aCenter_f /= sum_weights
        aCenter_m /= sum_weights
    # Subtract the centers-of-mass from the original coordinates for each object
    aaXf = aaXf_orig - aCenter_f
    aaXm = aaXm_orig - aCenter_m

    # Calculate the "M" array from the Diamond paper (equation 16)
    M = np.matmul(aaXm.T, (aaXf * aWeights))
    """ # For reference, here's the same code using for-loops
    M = np.zeros((3,3))
    for n in range(0, N):
        for i in range(0, 3):
            for j in range(0, 3):
                M[i,j] += aWeights[n] * aaXm[n,i] * aaXf[n,j]
    """

    # Calculate Q (equation 17)
    Q = M + M.T - 2*np.eye(3)*np.trace(M)

    # Calculate V (equation 18)
    V = np.empty(3)
    V[0] = M[1,2] - M[2,1]
    V[1] = M[2,0] - M[0,2]
    V[2] = M[0,1] - M[1,0]

    # Calculate "P" (equation 22)
    P = np.zeros((4,4))
    P[:3, :3] = Q  # copy Q into the first 3 rows and first 3 columns of P
    P[3, :3] = V   # copy V into the 4th row of P
    P[:3, 3] = V   # copy V into the 4th column of P


    # Calculate "p".  
    # "p" contains the optimal rotation (in backwards-quaternion format)
    # (Note: A discussion of various quaternion conventions is included below.)
    # First, specify the default value for p:
    p = np.zeros(4)
    p[3] = 1.0           # p = [0,0,0,1]    default value
    pPp: float = 0.0            # = p^T * P * p    (zero by default)
    singular: bool = (N < 2)   # (it doesn't make sense to rotate a single point)

    try:
        #http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.eigh.html
        aEigenvals, aaEigenvects = LA.eigh(P)

    except LA.LinAlgError:
        singular = True

    if not singular:  # (don't crash if the caller supplies nonsensical input)
        i_eval_max = np.argmax(aEigenvals)
        pPp = np.max(aEigenvals)
        p[:] = aaEigenvects[:, i_eval_max]
        """ # For reference, here's the same code using for-loops
        eval_max = aEigenvals[0]
        i_eval_max = 0
        for i in range(1, 4):
            if aEigenvals[i] > eval_max:
                eval_max = aEigenvals[i]
                i_eval_max = i
        p[0] = aaEigenvects[0,i_eval_max]
        p[1] = aaEigenvects[1,i_eval_max]
        p[2] = aaEigenvects[2,i_eval_max]
        p[3] = aaEigenvects[3,i_eval_max]
        pPp = eval_max
        """

    # normalize the vector
    # (It should be normalized already, but just in case it is not, do it again)
    p /= np.linalg.norm(p)

    # Finally, calculate the rotation matrix corresponding to "p"
    # (p is in backwards-quaternion format)

    aaRotate = np.empty((3,3))
    aaRotate[0,0] =  (p[0]*p[0])-(p[1]*p[1])-(p[2]*p[2])+(p[3]*p[3])
    aaRotate[1,1] = -(p[0]*p[0])+(p[1]*p[1])-(p[2]*p[2])+(p[3]*p[3])
    aaRotate[2,2] = -(p[0]*p[0])-(p[1]*p[1])+(p[2]*p[2])+(p[3]*p[3])
    aaRotate[0,1] = 2*(p[0]*p[1] - p[2]*p[3])
    aaRotate[1,0] = 2*(p[0]*p[1] + p[2]*p[3])
    aaRotate[1,2] = 2*(p[1]*p[2] - p[0]*p[3])
    aaRotate[2,1] = 2*(p[1]*p[2] + p[0]*p[3])
    aaRotate[0,2] = 2*(p[0]*p[2] + p[1]*p[3])
    aaRotate[2,0] = 2*(p[0]*p[2] - p[1]*p[3])

    ## Alternatively, in modern python versions, this code also works:
    # from scipy.spatial.transform import Rotation as R
    # the_rotation = R.from_quat(p)
    # aaRotate = the_rotation.as_matrix()

    # Optional: Decide the scale factor, c
    c: float = 1.0   # by default, don't rescale the coordinates
    if allow_rescale and (not singular):
        Waxaixai: float = np.sum(aWeights * aaXm * aaXm)
        WaxaiXai: float = np.sum(aWeights * aaXm * aaXf)
        """ # For reference, here's the same code using for-loops
        Waxaixai = 0.0
        WaxaiXai = 0.0
        for a in range(0, N):
            for i in range(0, 3):
                Waxaixai += aWeights[a,0] * aaXm[a,i] * aaXm[a,i]
                WaxaiXai += aWeights[a,0] * aaXm[a,i] * aaXf[a,i]
        """
        c = (WaxaiXai + pPp) / Waxaixai

    # Finally compute the RMSD between the two coordinate sets:
    # First compute E0 from equation 24 of the paper

    E0: float = np.sum(aWeights * (aaXf - c*aaXm)**2)
    sum_sqr_dist: float = max(0, E0 - c * 2.0 * pPp)
    """ # For reference, here is the same code using for-loops
    E0 = 0.0
    for n in range(0, N):
        for d in range(0, 3):
            E0 += aWeights[n] * ((aaXf[n,d] - c*aaXm[n,d])**2)
    sum_sqr_dist = E0 - c*2.0*pPp
    if sum_sqr_dist < 0.0: #(edge case due to rounding error)
        sum_sqr_dist = 0.0
    """

    rmsd: float = 0.0
    if sum_weights != 0.0:
        rmsd = np.sqrt(sum_sqr_dist/sum_weights)

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

    # new code (avoiding for-loops)
    aTranslate: NDArray = aCenter_f - np.matmul(c*aaRotate,aCenter_m).T.reshape(3,)
    """ # For reference, here is the same code using for-loops
    aTranslate = np.empty(3)
    for i in range(0,3):
        aTranslate[i] = aCenter_f[i]
        for j in range(0,3):
            aTranslate[i] -= c*aaRotate[i,j]*aCenter_m[j]
    """

    if report_quaternion: # does the caller want the quaternion?
        # The p array is a quaternion that uses this convention:
        #https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.transform.Rotation.from_quat.html
        # However it seems that the following convention is much more popular:
        # https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
        # https://mathworld.wolfram.com/Quaternion.html
        # So I return "q" (a version of "p" using the more popular convention).
        q = np.empty(4)
        q[0] = p[3]
        q[1] = p[0]
        q[2] = p[1]
        q[3] = p[2]
        return rmsd, q, aTranslate, c
    else:
        return rmsd, aaRotate, aTranslate, c
