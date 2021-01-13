#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2015 Huy Nguyen <huy.nguyendinh09@gmail.com>
#
# This file is part of python-cope.
#
# python-cope is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# python-cope is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# python-cope. If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import scipy.linalg
# Plots
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def TransformInv(T):
  """
  Calculates the inverse of the input homogeneous transformation.
  
  This method is more efficient than using C{numpy.linalg.inv}, given 
  the special properties of the homogeneous transformations.
  
  @type T: array, shape (4,4)
  @param T: The input homogeneous transformation
  @rtype: array, shape (4,4)
  @return: The inverse of the input homogeneous transformation
  """
  R = T[:3,:3].T
  p = T[:3,3]
  T_inv = np.identity(4)
  T_inv[:3,:3] = R
  T_inv[:3,3] = np.dot(-R, p)
  return T_inv

def TranValidate(T):
  """
  Validate T
  @type T:    array 4x4 
  @param T:   transformation matrix
  """
  raise NotImplementedError


def RotValidate(C):
  raise NotImplementedError


def TranAd(T):
  """
  Compute Adjoint of 4x4 transformation matrix, return a 6x6 matrix
  @type T:    array 4x4 
  @param T:   transformation matrix
  """
  C = T[:3,:3]
  r = T[:3,3]
  AdT = np.zeros([6,6])
  AdT[:3,:3] = C
  AdT[:3,3:] = np.dot(Hat(r),C)
  AdT[3:,3:] = C
  return AdT


def Hat(vec):
  """
  hat operator - return skew matrix (return 3x3 or 4x4 matrix) from vec
  @param vec:   vector of 3 (rotation) or 6 (transformation)
  """
  if vec.shape[0] == 3: # skew from vec
    return np.array([[0,-vec[2],vec[1]],[vec[2],0,-vec[0]],[-vec[1],vec[0],0]])
  elif vec.shape[0] == 6:
    vechat = np.zeros((4,4))
    vechat[:3,:3] = Hat(vec[3:])
    vechat[:3,3] = vec[:3]
    return vechat
  else:
    raise ValueError("Invalid vector length for hat operator\n")


def VecFromSkew(r):
  return np.array([r[2,1],r[0,2],r[1,0]])


def CurlyHat(vec):
  """
  Builds the 6x6 curly hat matrix from the 6x1 input
  @param vec:          a 6x1 vector xi
  @param veccurlyhat:  a 6x6 matrix 
  """
  veccurlyhat = np.zeros((6,6))
  veccurlyhat[:3,:3] = Hat(vec[3:])
  veccurlyhat[:3,3:] = Hat(vec[:3])
  veccurlyhat[3:,3:] = Hat(vec[3:])
  return veccurlyhat


def CovOp1(A):
  """ 
  Covariance operator 1 - eq. 44
  """
  return -np.trace(A)*np.eye(3) + A

 
def CovOp2(A,B):
  """ 
  Covariance operator 2 - eq. 45
  """
  return np.dot(CovOp1(A),CovOp1(B)) + CovOp1(np.dot(B,A))


def TranToVec(T):
  """
  Compute the matrix log of the transformation matrix T
  Convert from T to xi
  @param T:       4x4
  @param return:  return a 6x1 vector in tangent coordinates computed from T.
  """
  C = T[:3,:3]
  r = T[:3,3]
  
  phi = RotToVec(C)
  invJ = VecToJacInv(phi)
  
  rho = np.dot(invJ,r)
  return np.hstack([rho,phi])


# def RotToVec(C):
#   """
#   Compute the matrix log of the rotation matrix C
#   @param C:      3x3
#   @param return: Return a 3x1 vector (axis*angle) computed from C
#   """
#   #rotValidate(C)
#   if(abs(np.trace(C)+1)>1e-10):
#     if(np.linalg.norm(C-np.eye(3))<=1e-10):
#       return np.zeros(3)
#     else:
#       phi = np.arccos((np.trace(C)-1)/2)
#       return VecFromSkew(phi/(2*np.sin(phi))*(C-C.T))
#   else:
#     eigval, eigvect = np.linalg.eig(C)
#     for (i,val) in enumerate(eigval):
#       if abs((val-1)) <= 1e-10:
#         return np.pi*np.real(eigvect[:,i])

def RotToVec(C):
  """
  Compute the matrix log of the rotation matrix C
  @param C:      3x3
  @param return: Return a 3x1 vector (axis*angle) computed from C
  """
  # RotValidate(C)
  epsilon = 0.0001
  epsilon2 = 0.001
  if ((abs(C[0,1]-C[1,0])<epsilon) and (abs(C[0,2]-C[2,0])<epsilon) and (abs(C[1,2]-C[2,1])<epsilon)):
    # singularity found
    # first check for identity matrix which must have +1 for all terms
		# in leading diagonaland zero in other terms
    if ((abs(C[0,1]+C[1,0]) < epsilon2) and (abs(C[0,2]+C[2,0]) < epsilon2) and (abs(C[1,2]+C[2,1]) < epsilon2) and (abs(C[0,0]+C[1,1]+C[2,2]-3) < epsilon2)): # this singularity is identity matrix so angle = 0
      return np.zeros(3) #zero angle, arbitrary axis 
    # otherwise this singularity is angle = 180
    angle = np.pi
    xx = (C[0,0]+1)/2.
    yy = (C[1,1]+1)/2.
    zz = (C[2,2]+1)/2.
    xy = (C[0,1]+C[1,0])/4.
    xz = (C[0,2]+C[2,0])/4.
    yz = (C[1,2]+C[2,1])/4.
    if ((xx > yy) and (xx > zz)): # C[0][0] is the largest diagonal term
      if (xx< epsilon):
        x = 0
        y = np.sqrt(2)/2.
        z = np.sqrt(2)/2.
      else:
        x = np.sqrt(xx)
        y = xy/x
        z = xz/x
    elif (yy > zz): # C[1][1] is the largest diagonal term
      if (yy< epsilon):
        x = np.sqrt(2)/2.
        y = 0
        z = np.sqrt(2)/2.
      else:
        y = np.sqrt(yy)
        x = xy/y
        z = yz/y
    else: # C[2][2] is the largest diagonal term so base result on this
      if (zz< epsilon):
        x = np.sqrt(2)/2.
        y = np.sqrt(2)/2.
        z = 0
      else:
        z = np.sqrt(zz)
        x = xz/z
        y = yz/z
    return angle*np.array((x,y,z))
  s = np.sqrt((C[2,1] - C[1,2])*(C[2,1] - C[1,2])+(C[0,2] - C[2,0])*(C[0,2] - C[2,0])+(C[1,0] - C[0,1])*(C[1,0] - C[0,1])) # used to normalise
  if (abs(s) < 0.001):
    # prevent divide by zero, should not happen if matrix is orthogonal and should be
    # caught by singularity test above, but I've left it in just in case
    s=1 
        
  angle = np.arccos(( C[0,0] + C[1,1] + C[2,2] - 1)/2.)
  x = (C[2,1] - C[1,2])/s
  y = (C[0,2] - C[2,0])/s
  z = (C[1,0] - C[0,1])/s
  return angle*np.array((x,y,z))

def VecToRot(phi):
  """
  Return a rotation matrix computed from the input vec (phi 3x1)
  @param phi: 3x1 vector (input)
  @param C:   3x3 rotation matrix (output)
  """
  tiny = 1e-12
  #check for small angle
  nr = np.linalg.norm(phi)
  if nr < tiny:
    #~ # If the angle (nr) is small, fall back on the series representation.
    # C = VecToRotSeries(phi,10)
    C = np.eye(3)
  else:
    R = Hat(phi)
    C = np.eye(3) + np.sin(nr)/nr*R + (1-np.cos(nr))/(nr*nr)*np.dot(R,R)
  return C


def VecToRotSeries(phi, N):
  """"
  Build a rotation matrix using the exponential map series with N elements in the series 
  @param phi: 3x1 vector
  @param N:   number of terms to include in the series
  @param C:   3x3 rotation matrix (output)
  """
  C = np.eye(3)
  xM = np.eye(3)
  cmPhi = Hat(phi)
  for n in range(N):
    xM = np.dot(xM, cmPhi)/(n+1)
    C = C + xM
  # Project the resulting rotation matrix back onto SO(3)
  C = np.dot(C,np.linalg.inv(scipy.linalg.sqrtm(np.dot(C.T,C))))
  return C


def cot(x):
  return 1./np.tan(x)

  
def VecToJacInv(vec):
  """
  Construction of the 3x3 J^-1 matrix or 6x6 J^-1 matrix.
  @param vec:  3x1 vector or 6x1 vector (input)
  @param invJ: 3x3 inv(J) matrix or 6x6 inv(J) matrix (output)
  """
  tiny = 1e-12
  if vec.shape[0] == 3: # invJacobian of SO3
    phi = vec
    nr = np.linalg.norm(phi)
    if nr < tiny:
      # If the angle is small, fall back on the series representation
      invJSO3 = VecToJacInvSeries(phi,10)
    else:
      axis = phi/nr
      invJSO3 = 0.5*nr*cot(nr*0.5)*np.eye(3) + (1- 0.5*nr*cot(0.5*nr))*axis[np.newaxis]*axis[np.newaxis].T- 0.5*nr*Hat(axis)
    return invJSO3
  elif vec.shape[0] == 6: # invJacobian of SE3
    rho = vec[:3]
    phi = vec[3:]
    
    nr = np.linalg.norm(phi)
    if nr < tiny:
      # If the angle is small, fall back on the series representation
      invJSO3 = VecToJacInvSeries(phi,10)
    else:
      invJSO3 = VecToJacInv(phi)
    Q = VecToQ(vec)
    invJSE3 = np.zeros((6,6))
    invJSE3[:3,:3] = invJSO3
    invJSE3[:3,3:] = -np.dot(np.dot(invJSO3,Q), invJSO3)
    invJSE3[3:,3:] = invJSO3
    return invJSE3
  else:
    raise ValueError("Invalid input vector length\n")


def VecToJacInvSeries(vec,N):
  """
  Construction of the 3x3 J^-1 matrix or 6x6 J^-1 matrix. Series representation.
  @param vec:  3x1 vector or 6x1 vector
  @param N:    number of terms to include in the series
  @param invJ: 3x3 inv(J) matrix or 6x6 inv(J) matrix (output)
  """
  if vec.shape[0] == 3: # invJacobian of SO3
    invJSO3 = np.eye(3)
    pxn = np.eye(3)
    px = Hat(vec)
    for n in range(N):
      pxn = np.dot(pxn,px)/(n+1)
      invJSO3 = invJSO3 + BernoulliNumber(n+1)*pxn
    return invJSO3
  elif vec.shape[0] == 6: # invJacobian of SE3
    invJSE3 =np.eye(6)
    pxn = np.eye(6)
    px = CurlyHat(vec)
    for n in range(N):
      pxn = np.dot(pxn,px)/(n+1)
      invJSE3 = invJSE3 + BernoulliNumber(n+1)*pxn
    return invJSE3
  else:
    raise ValueError("Invalid input vector length\n")


def BernoulliNumber(n):
  """
  Generate Bernoulli number
  @param n:  interger (0,1,2,...)
  """
  from fractions import Fraction as Fr
  if n == 1: return -0.5
  A = [0] * (n+1)
  for m in range(n+1):
    A[m] = Fr(1, m+1)
    for j in range(m, 0, -1):
      A[j-1] = j*(A[j-1] - A[j])
  return A[0].numerator*1./A[0].denominator # (which is Bn)


def VecToJac(vec):
  """ 
  Construction of the J matrix
  @param vec: a 3x1 vector for SO3 or a 6x1 vector for SE3 (input)
  @param J:   a 3x3 J matrix for SO3 or a 6x6 J matrix for SE3 (output)
  """
  tiny = 1e-12
  if vec.shape[0] == 3: # Jacobian of SO3
    phi = vec
    nr = np.linalg.norm(phi)
    if nr < tiny:
      # If the angle is small, fall back on the series representation
      JSO3 = VecToJacSeries(phi,10)
    else:
      axis = phi/nr
      cnr = np.cos(nr)
      snr = np.sin(nr)
      JSO3 = (snr/nr)*np.eye(3) + (1-snr/nr)*axis[np.newaxis]*axis[np.newaxis].T + ((1-cnr)/nr)*Hat(axis)
    return JSO3
  elif vec.shape[0] == 6: # Jacobian of SE3
    rho = vec[:3]
    phi = vec[3:]
    nr = np.linalg.norm(phi)
    if nr < tiny:
      #If the angle is small, fall back on the series representation
      JSE3 = VecToJacSeries(phi,10);
    else:
      JSO3 = VecToJac(phi)
      Q = VecToQ(vec)
      JSE3 = np.zeros((6,6))
      JSE3[:3,:3] = JSO3
      JSE3[:3,3:] = Q
      JSE3[3:,3:] = JSO3
    return JSE3
  else:
    raise ValueError("Invalid input vector length\n")


def VecToJacSeries(vec,N):
  """ 
  Construction of the J matrix from Taylor series
  @param vec: a 3x1 vector for SO3 or a 6x1 vector for SE3 (input)
  @param N:   number of terms to include in the series (input)
  @param J:   a 3x3 J matrix for SO3 or a 6x6 J matrix for SE3 (output)
  """
  if vec.shape[0] == 3: # Jacobian of SO3
    JSO3 = np.eye(3)
    pxn = np.eye(3)
    px = Hat(vec)
    for n in range(N):
      pxn = np.dot(pxn,px)/(n+2)
      JSO3 = JSO3 + pxn
    return JSO3
  elif vec.shape[0] == 6: # Jacobian of SE3
    JSE3 = np.eye(6)
    pxn = np.eye(6)
    px = CurlyHat(vec)
    for n in range(N):
      pxn = np.dot(pxn,px)/(n+2)
      JSE3 = JSE3 + pxn
    return JSE3
  else:
    raise ValueError("Invalid input vector length\n")
  return


def VecToQ(vec):
  """
  @param vec: a 6x1 vector (input)
  @param Q:   the 3x3 Q matrix (output)
  """
  rho = vec[:3]
  phi = vec[3:]
  
  nr = np.linalg.norm(phi)
  if nr == 0:
    nr = 1e-12
  nr2 = nr*nr
  nr3 = nr2*nr
  nr4 = nr3*nr
  nr5 = nr4*nr
  
  cnr = np.cos(nr)
  snr = np.sin(nr)
  
  rx = Hat(rho)
  px = Hat(phi)
  
  t1 = 0.5*rx
  t2 = ((nr-snr)/nr3)*(np.dot(px,rx) + np.dot(rx,px) + np.dot(np.dot(px,rx),px))
  m3 = (1-0.5*nr2 - cnr)/nr4
  t3 = -m3*(np.dot(np.dot(px,px),rx) +np.dot(np.dot(rx,px),px) -3*np.dot(np.dot(px,rx),px))
  m4 = 0.5*(m3 - 3*(nr - snr -nr3/6)/nr5)
  t4 = -m4*(np.dot(px,np.dot(np.dot(rx,px),px)) + np.dot(px,np.dot(np.dot(px,rx),px)))
  Q = t1 + t2 + t3 + t4
  return Q


def VecToTran(vec):
  """
  Build a transformation matrix using the exponential map, closed form
  @param vec: 6x1 vector (input)
  @param T:   4x4 transformation matrix (output)
  """
  rho = vec[:3]
  phi = vec[3:]
  
  C = VecToRot(phi)
  JSO3 = VecToJac(phi)
  
  T = np.eye(4)
  T[:3,:3] = C
  T[:3,3] = np.dot(JSO3,rho)
  return T


def VecToTranSeries(p, N):
  """
  Build a transformation matrix using the exponential map series with N elements in the series
  @param p: 6x1 vector (input)
  @param N: number of terms to include in the series (input)
  @param T: 4x4 transformation matrix (output)
  """
  T = np.eye(4)
  xM = np.eye(4)
  bpP = Hat(p)
  for n in range(N):
    xM = np.dot(xM, bpP/(n+1))
    T = T + xM
  return T


def Propagating(T1, sigma1, T2, sigma2, method = 2):
  """
  Find the total uncertainty in a compound spatial relation (Compounding two uncertain transformations)
  @param T1:     4x4 mean of left transformation 
  @param sigma1: 6x6 covariance of left transformation
  @param T2:     4x4 mean of right transformation
  @param sigma2: 6x6 covariance of right transformations
  @param method: an integer indicating the method to be used to perform compounding
                 (1=second-order, 2=fourth-order)
  @param T:      4x4 mean of compounded transformation (output)
  @param sigma:  6x6 covariance of compounded transformation (output)
  """
  # Compound the means
  T = np.dot(T1,T2)
  # Compute Adjoint of transformation T1
  AdT1 = TranAd(T1)
  sigma2prime = np.dot(np.dot(AdT1,sigma2),AdT1)
  if method == 1:
    # Second-order method
    sigma = sigma1 + sigma2prime    
  elif method == 2:
    # Fourth-order method
    sigma1rr = sigma1[:3,:3]
    sigma1rp = sigma1[:3,3:]
    sigma1pp = sigma1[3:,3:]
    
    sigma2rr = sigma2prime[:3,:3]
    sigma2rp = sigma2prime[:3,3:]
    sigma2pp = sigma2prime[3:,3:]
    
    A1 = np.zeros((6,6))
    A1[:3,:3] = CovOp1(sigma1pp)
    A1[:3,3:] = CovOp1(sigma1rp + sigma1rp.T)
    A1[3:,3:] = CovOp1(sigma1pp)
    
    A2 = np.zeros((6,6))
    A2[:3,:3] = CovOp1(sigma2pp)
    A2[:3,3:] = CovOp1(sigma2rp + sigma2rp.T)
    A2[3:,3:] = CovOp1(sigma2pp)

    Brr = CovOp2(sigma1pp,sigma2rr) + CovOp2(sigma1rp.T,sigma2rp) + CovOp2(sigma1rp,sigma2rp.T) + CovOp2(sigma1rr,sigma2pp)
    Brp = CovOp2(sigma1pp,sigma2rp.T) + CovOp2(sigma1rp.T,sigma2pp)
    Bpp = CovOp2(sigma1pp, sigma2pp)
    
    B = np.zeros((6,6))
    B[:3,:3] = Brr
    B[:3,3:] = Brp
    B[3:,:3] = Brp.T
    B[3:,3:] = Bpp
    
    sigma = sigma1 + sigma2prime + 1/12.*(np.dot(A1,sigma2prime)+np.dot(sigma2prime,A1.T) + np.dot(sigma1,A2) + np.dot(sigma1,A2.T)) + B/4.   
  return T, sigma

def PropagatingWithSeparateRotTrans(R1,sigmaR1,t1,sigmat1,R2,sigmaR2,t2,sigmat2):
    """
    Find the total uncertainty in a compound spatial relation (Compounding two uncertain transformations) where we separate rotation and translation
    @param sigmaR1,sigmaR2: cov of Rot vec (RotToVec(T[:3,:3]))
    @param sigmat1,sigamt2: cov of Trans vec (T[:3,3])
    """
    # Compound the means
    R = np.dot(R1,R2)
    t = np.dot(R1,t2)+t1
    # Compute the cov of the compounding rot
    sigmaR2prime = np.dot(np.dot(R1,sigmaR2),np.transpose(R1))
    A1 = CovOp1(sigmaR1)
    A2 = CovOp1(sigmaR2prime)
    B = CovOp2(sigmaR1, sigmaR2prime)
    sigmaR = sigmaR1 + sigmaR2prime + 1/12.*(np.dot(A1,sigmaR2prime) + np.dot(sigmaR2prime,A1.T) + np.dot(sigmaR1,A2) + np.dot(sigmaR1,A2.T)) + B/4.
    # Compute the cov of teh compounding trans
    R1t2 = np.dot(R1,t2)
    sigmat = sigmat1 + np.dot(np.dot(R1,sigmat2),np.transpose(R1)) + np.dot(np.dot(Hat(R1t2),sigmaR1),np.transpose(Hat(R1t2)))
    return R, sigmaR, t, sigmat


def Fusing(Tlist, sigmalist, N = 0, maxiterations=30, retiter=False):
  """
  Find the total uncertainty in a compound spatial relation (Compounding two uncertain transformations)
  @param Tlist:     a list of 4x4 transformations
  @param sigmalist: a list of corresponding 6x6 covariance matrices
  @param N:         N == 0(default):JacInv is computed analytically using eq. 100
                    N != 0: JacInv is computed using eq. 103, using N first terms in the eq.
  @param T:      4x4 mean of fused transformation (output)
  @param sigma:  6x6 covariance of fused transformation (output)
  """
  assert len(Tlist) == len(sigmalist), "Invalid data list length\n"
  kmax = len(Tlist)
  
  T = Tlist[0]
  Vprv = 0
  for i in range(maxiterations): # Gauss-Newton iterations
    LHS = np.zeros(6)
    RHS = np.zeros(6)
    for k in range(kmax):
      xik = TranToVec(np.dot(T,np.linalg.inv(Tlist[k])))
      if N ==0:
        invJ = VecToJacInv(xik)
      else:
        invJ = VecToJacInvSeries(xik, N)
      invJtS = np.dot(invJ.T, np.linalg.inv(sigmalist[k]))
      LHS = LHS + np.dot(invJtS,invJ)
      RHS = RHS + np.dot(invJtS, xik)
    xi = -np.linalg.solve(LHS,RHS)
    # print "xi", xi
    T = np.dot(VecToTran(xi),T)
    # print "T", T
    sigma = np.linalg.inv(LHS)
    # How low did the objective function get?
    V = 0.
    for k in range(kmax):
      xik = TranToVec(np.dot(T,np.linalg.inv(Tlist[k])))
      V = V + np.dot(np.dot(xik.T,np.linalg.inv(sigmalist[k])),xik) / 2.
    if abs(V - Vprv) < 1e-10:
      break 
    Vprv = V
  if retiter:
    return T, sigma, i+1
  else:
    return T, sigma

def CovInverseTran(T,sigma):
    """
    Compute the cov of the inverse transformation. (Follow Ethan Eade's note on lie group.)
    """
    Tinv = np.linalg.inv(T)
    AdTinv = TranAd(Tinv)
    sigmaTinv = np.dot(np.dot(AdTinv,sigma),np.transpose(AdTinv))
    return Tinv, sigmaTinv

def CovInverseTranWithSeparateRotTrans(R,sigmaR,t,sigmat):
    """
    Compute the cov of the inverse transformation where Rot and Trans 's noises are assumed to be independent
    """
    Rinv = np.linalg.inv(R)
    tinv = -np.dot(Rinv,t)
    sigmaRinv = sigmaR
    hatRinvt = Hat(np.dot(Rinv,t))
    sigmatinv = np.dot(np.dot(hatRinvt,sigmaRinv),np.transpose(hatRinvt)) + np.dot(np.dot(Rinv,sigmat),R)
    return Rinv, sigmaRinv, tinv, sigmatinv

def Visualize(Tlist,sigmalist, nsamples = 100):
  """
  Visualize an estimation (a point will be used to represent the translation position of a transformation)
  @param Tlist:     a list of Transformations
  @param sigmalist: a list of corresponding sigmas
  @param nsamples:  the number of samples generated for each (T,sigma)
  """
  import matplotlib.cm as cm
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  cholsigmalist = []
  colors = iter(cm.rainbow(np.linspace(0, 1, len(Tlist))))
  for i in range(len(sigmalist)):
    color = next(colors)
    cholsigma = np.linalg.cholesky(sigmalist[i]).T
    Tsample = []
    for k in range(nsamples):
      vecsample = np.dot(cholsigma,np.random.randn(6,1))
      #vecsample = np.dot(cholsigma, np.random.uniform(-1,1,size = 6))
      vecsample.resize(6)
      Tsample = np.dot(VecToTran(vecsample), Tlist[i])
      ax.scatter(Tsample[0,3],Tsample[1,3],Tsample[2,3], c = color)

  ax.set_autoscaley_on(False)
  ax.set_xlim([-0.5, 0.5])
  ax.set_ylim([-0.5, 0.5])
  ax.set_zlim([-0.5, 0.5])
  ax.set_xlabel('X Label')
  ax.set_ylabel('Y Label')
  ax.set_zlabel('Z Label')
  plt.show(False)
  return True


def IsInside(point, center, sigma):
  """
  Check whether a point (transformation) is in the region formed by (center,sigma) or not 
  """
  cholsigma = np.linalg.cholesky(sigma).T
  univariable = np.dot(np.linalg.inv(cholsigma),(point-center))
  nr = np.linalg.norm(univariable)
  if nr <= 1.0:
    return True
  else:
    return False

class Pose():
  # SE3 pose where Rot and Trans are separated.
  def __init__(self,rot, sigmarot, trans,sigmatrans):
    self.rot = rot
    self.trans = trans
    self.sigmarot = sigmarot
    self.sigmatrans = sigmatrans
    self.transform = np.eye(4)
    self.transform[:3,:3] = rot
    self.transform[:3,3] = trans

  # def transform(self):
  #   transform = np.eye(4)
  #   transform[:3,:3] = self.rot
  #   transform[:3,3] = self.trans
  #   return transform

def ConstPose(T):
  '''
  Return a pose obj w no covariance matrix. input: a 4x4 matrix
  '''
  sigmarot = np.zeros((3,3))
  sigmatrans = np.zeros((3,3))
  return Pose(T[:3,:3], sigmarot, T[:3,3],sigmatrans)

def UpdatePose(pose):
  pose.rot = pose.transform[:3,:3]
  pose.trans = pose.transform[:3,3]
  return True

def Dot(pose1,pose2):
  '''Find the total uncertainty in a compound spatial relation (Compounding two uncertain transformations) where we separate rotation and translation.
  output: a pose Pose(R, sigmaR, t, sigmat)
  '''
  UpdatePose(pose1)
  UpdatePose(pose2)
  R, sigmaR, t, sigmat = PropagatingWithSeparateRotTrans(pose1.rot,pose1.sigmarot,pose1.trans,pose1.sigmatrans,pose2.rot,pose2.sigmarot,pose2.trans,pose2.sigmatrans)
  return Pose(R, sigmaR, t, sigmat)

def Inverse(pose):
  '''
  Return the inverse and cov of the the inverse transformation
  '''
  UpdatePose(pose)
  R, sigmaR, t, sigmat = CovInverseTranWithSeparateRotTrans(pose.rot,pose.sigmarot,pose.trans,pose.sigmatrans)
  return Pose(R, sigmaR, t, sigmat)
