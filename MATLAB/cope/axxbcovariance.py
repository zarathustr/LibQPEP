#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2017 Huy Nguyen <huy.nguyendinh09@gmail.com>
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

import copy
import math
import numpy as np
import random

import cope.SE3lib as SE3

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
plt.ion()

def Eigsorted(cov):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]

def VisualizeCovariances(cov_rot, cov_trans, minx,maxx,miny,maxy):
  plt.subplot(231) # x and y axis
  nstd=1
  alpha = 0.5
  mean = (0,0)
  cov0 = cov_rot

  cov = [[cov0 [0][0],cov0 [0][1]],[cov0 [1][0],cov0 [1][1]]]
  pos = mean

  ax = plt.gca()
  vals, vecs = Eigsorted(cov)
  theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
  # Width and height are "full" widths, not radius
  width, height = 2 * nstd * np.sqrt(vals)
  ellip1 = Ellipse(xy=pos, width=width, height=height, angle=theta,alpha=0.5,color='red', linewidth=2, fill=False)
  ax.add_artist(ellip1)
  plt.axis([minx,maxx,miny,maxy])
  plt.xlabel(r'${\bf{\xi}}_{\bf{R} x} (rad)$',fontsize=20, labelpad=8)
  plt.ylabel(r'${\bf{\xi}}_{\bf{R} y} (rad)$',fontsize=20, labelpad=-8)
  plt.xticks(np.arange(minx, maxx+maxx/2, (maxx-minx)/2))
  plt.yticks(np.arange(miny, maxy+maxy/2, (maxy-miny)/2))
  ax.set(aspect='equal')


  plt.subplot(234)
  mean = (0,0)
  cov = cov_trans
  cov = [[cov [0][0],cov [0][1]],[cov [1][0],cov [1][1]]]
  pos=mean
  ax = plt.gca()
  vals, vecs = Eigsorted(cov)
  theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
  # Width and height are "full" widths, not radius
  width, height = 2 * nstd * np.sqrt(vals)
  ellip2 = Ellipse(xy=pos, width=width, height=height, angle=theta,alpha=0.5,color='red', linewidth=2, fill=False)
  ax.add_artist(ellip2)
  plt.axis([minx,maxx,miny,maxy])
  plt.xlabel(r'${\bf{\xi}}_{\bf{t} x} (mm)$',fontsize=20, labelpad=8)
  plt.ylabel(r'${\bf{\xi}}_{\bf{t} y} (mm)$',fontsize=20, labelpad=-8)
  plt.xticks(np.arange(minx, maxx+maxx/2, (maxx-minx)/2))
  plt.yticks(np.arange(miny, maxy+maxy/2, (maxy-miny)/2))
  ax.set(aspect='equal')
 ############################################################

  plt.subplot(232)
  cov = cov_rot
  cov = [[cov0 [1][1],cov0 [1][2]],[cov0 [2][1],cov0 [2][2]]]
  pos = mean

  ax = plt.gca()
  vals, vecs = Eigsorted(cov)
  theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
  # Width and height are "full" widths, not radius
  width, height = 2 * nstd * np.sqrt(vals)
  ellip1 = Ellipse(xy=pos, width=width, height=height, angle=theta,alpha=0.5,color='red', linewidth=2, fill=False)
  ax.add_artist(ellip1)
  plt.axis([minx,maxx,miny,maxy])
  plt.xlabel(r'${\bf{\xi}}_{\bf{R} y} (rad)$',fontsize=20, labelpad=8)
  plt.ylabel(r'${\bf{\xi}}_{\bf{R} z} (rad)$',fontsize=20, labelpad=-8)  
  plt.xticks(np.arange(minx, maxx+maxx/2, (maxx-minx)/2))
  plt.yticks(np.arange(miny, maxy+maxy/2, (maxy-miny)/2))
  ax.set(aspect='equal')

  plt.subplot(235)
  mean = (0,0)
  cov = cov_trans
  cov = [[cov [1][1],cov [1][2]],[cov [2][1],cov [2][2]]]
  pos=mean
  ax = plt.gca()
  vals, vecs = Eigsorted(cov)
  theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
  # Width and height are "full" widths, not radius
  width, height = 2 * nstd * np.sqrt(vals)
  ellip2 = Ellipse(xy=pos, width=width, height=height, angle=theta,alpha=0.5,color='red', linewidth=2, fill=False)
  ax.add_artist(ellip2)
  plt.axis([minx,maxx,miny,maxy])
  plt.xlabel(r'${\bf{\xi}}_{\bf{t} y} (mm)$',fontsize=20, labelpad=8)
  plt.ylabel(r'${\bf{\xi}}_{\bf{t} z} (mm)$',fontsize=20, labelpad=-8)
  plt.xticks(np.arange(minx, maxx+maxx/2, (maxx-minx)/2))
  plt.yticks(np.arange(miny, maxy+maxy/2, (maxy-miny)/2))
  ax.set(aspect='equal')

  ###################################################################

  plt.subplot(233)
  cov = cov_rot
  cov = [[cov0 [0][0],cov0 [0][2]],[cov0 [2][0],cov0 [2][2]]]
  pos = mean

  ax = plt.gca()
  vals, vecs = Eigsorted(cov)
  theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
  # Width and height are "full" widths, not radius
  width, height = 2 * nstd * np.sqrt(vals)
  ellip1 = Ellipse(xy=pos, width=width, height=height, angle=theta,alpha=0.5,color='red', linewidth=2, fill=False)
  ax.add_artist(ellip1)
  plt.axis([minx,maxx,miny,maxy])
  plt.xticks(np.arange(minx, maxx+maxx/2, (maxx-minx)/2))
  plt.yticks(np.arange(miny, maxy+maxy/2, (maxy-miny)/2))
  plt.xlabel(r'${\bf{\xi}}_{\bf{R} x} (rad)$',fontsize=20, labelpad=8)
  plt.ylabel(r'${\bf{\xi}}_{\bf{R} z} (rad)$',fontsize=20, labelpad=-8)
  ax.set(aspect='equal')
  plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.8, hspace=None)

  plt.subplot(236)
  mean = (0,0)
  cov = cov_trans
  cov = [[cov [0][0],cov [0][2]],[cov [2][0],cov [2][2]]]

  pos=mean
  ax = plt.gca()
  vals, vecs = Eigsorted(cov)
  theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
  # Width and height are "full" widths, not radius
  width, height = 2 * nstd * np.sqrt(vals)
  ellip2 = Ellipse(xy=pos, width=width, height=height, angle=theta,alpha=0.5,color='red', linewidth=2, fill=False)
  ax.add_artist(ellip2)
  plt.axis([minx,maxx,miny,maxy])
  plt.xticks(np.arange(minx, maxx+maxx/2, (maxx-minx)/2))
  plt.yticks(np.arange(miny, maxy+maxy/2, (maxy-miny)/2))
  plt.xlabel(r'${\bf{\xi}}_{\bf{t} x}(mm)$',fontsize=20, labelpad=8)
  plt.ylabel(r'${\bf{\xi}}_{\bf{t} z}(mm)$',fontsize=20, labelpad=-8)
  ax.set(aspect='equal')
  plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.8, hspace=None)
  # plt.legend(handles=[ellip1],loc='upper center',ncol=2, bbox_to_anchor=(0.5,1.5))
  return True

def VisualizeRealEstCov(cov_real, cov_est, minx,maxx,miny,maxy,param):
    if param=='rot':
        subplotnum = 230
    if param=='trans':
        subplotnum = 233
    # Compare y z
    nstd=1
    alpha = 0.5
    mean = (0,0)
    cov0 = cov_real

    plt.subplot(subplotnum+1)
    cov = [[cov0 [0][0],cov0 [0][1]],[cov0 [1][0],cov0 [1][1]]]
    pos = mean

    ax = plt.gca()
    vals, vecs = Eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)
    ellip1 = Ellipse(xy=pos, width=width, height=height, angle=theta,alpha=0.5,color='green',linewidth=2,linestyle='dashed', fill=False,label='Empirical Estimation')
    ax.add_artist(ellip1)
    
    mean = (0,0)
    cov = cov_est
    cov = [[cov [0][0],cov [0][1]],[cov [1][0],cov [1][1]]]

    pos=mean
    ax = plt.gca()
    vals, vecs = Eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)
    ellip2 = Ellipse(xy=pos, width=width, height=height, angle=theta,alpha=0.5,color='red', linewidth=2, fill=False,label='Our algorithm')
    ax.add_artist(ellip2)
    ###########real data
    # plt.axis([-0.0025,0.0025,-0.002,0.0025])
    # plt.axis([-0.018,0.018,-0.018,0.018])
    ###########synthetic
    # plt.axis([-0.00055,0.00055,-0.00055,0.00055])
    plt.axis([minx,maxx,miny,maxy])
    if param=='rot':
        plt.xlabel(r'${\bf{\xi}}_{\bf{R} x} (rad)$',fontsize=20, labelpad=8)
        plt.ylabel(r'${\bf{\xi}}_{\bf{R} y} (rad)$',fontsize=20, labelpad=-8)
    if param=='trans':
        plt.xlabel(r'${\bf{\xi}}_{\bf{t} x} (mm)$',fontsize=20, labelpad=8)
        plt.ylabel(r'${\bf{\xi}}_{\bf{t} y} (mm)$',fontsize=20, labelpad=-8)
    plt.xticks(np.arange(minx, maxx+maxx/2, (maxx-minx)/2))
    plt.yticks(np.arange(miny, maxy+maxy/2, (maxy-miny)/2))
    ax.set(aspect='equal')

# ################################################################

    plt.subplot(subplotnum+2)
    cov = [[cov0 [1][1],cov0 [1][2]],[cov0 [2][1],cov0 [2][2]]]
    pos = mean

    ax = plt.gca()
    vals, vecs = Eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)
    ellip1 = Ellipse(xy=pos, width=width, height=height, angle=theta,alpha=0.5,color='green',linewidth=2,linestyle='dashed', fill=False,label='Empirical Estimation')
    ax.add_artist(ellip1)
    
    mean = (0,0)
    cov = cov_est
    cov = [[cov [1][1],cov [1][2]],[cov [2][1],cov [2][2]]]

    pos=mean
    ax = plt.gca()
    vals, vecs = Eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)
    ellip2 = Ellipse(xy=pos, width=width, height=height, angle=theta,alpha=0.5,color='red', linewidth=2, fill=False,label='Our algorithm')
    ax.add_artist(ellip2)
    ############real data
    # plt.axis([-0.0025,0.0025,-0.002,0.0025])
    #plt.axis([-0.018,0.018,-0.018,0.018])
    ###############synthetic
    # plt.axis([-0.00055,0.00055,-0.00055,0.00055])
    plt.axis([minx,maxx,miny,maxy])
    if param=='rot':
        plt.xlabel(r'${\bf{\xi}}_{\bf{R} y} (rad)$',fontsize=20, labelpad=8)
        plt.ylabel(r'${\bf{\xi}}_{\bf{R} z} (rad)$',fontsize=20, labelpad=-8)
    if param=='trans':
        plt.xlabel(r'${\bf{\xi}}_{\bf{t} y} (mm)$',fontsize=20, labelpad=8)
        plt.ylabel(r'${\bf{\xi}}_{\bf{t} z} (mm)$',fontsize=20, labelpad=-8)
    # plt.legend(handles=[ellip1, ellip2])
    plt.xticks(np.arange(minx, maxx+maxx/2, (maxx-minx)/2))
    plt.yticks(np.arange(miny, maxy+maxy/2, (maxy-miny)/2))
    if param == 'rot':
        plt.legend(handles=[ellip1, ellip2],loc='upper center',ncol=2, bbox_to_anchor=(0.5,1.5))
    ax.set(aspect='equal')

    ###################################################################

    plt.subplot(subplotnum+3)
    cov = [[cov0 [0][0],cov0 [0][2]],[cov0 [2][0],cov0 [2][2]]]
    pos = mean


    ax = plt.gca()
    vals, vecs = Eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)
    ellip1 = Ellipse(xy=pos, width=width, height=height, angle=theta,alpha=0.5,color='green',linewidth=2,linestyle='dashed', fill=False,label='Empirical Estimation')
    ax.add_artist(ellip1)

    mean = (0,0)
    cov = cov_est
    cov = [[cov [0][0],cov [0][2]],[cov [2][0],cov [2][2]]]

    pos=mean
    ax = plt.gca()
    vals, vecs = Eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)
    ellip2 = Ellipse(xy=pos, width=width, height=height, angle=theta,alpha=0.5,color='red', linewidth=2, fill=False,label ='Our algorithm')
    ax.add_artist(ellip2)
    ############real data
    # plt.axis([-0.0025,0.0025,-0.002,0.0025])
    # plt.axis([-0.018,0.018,-0.018,0.018])
    ###############synthetic
    # plt.axis([-0.00055,0.00055,-0.00055,0.00055])
    plt.axis([minx,maxx,miny,maxy])
    plt.xticks(np.arange(minx, maxx+maxx/2, (maxx-minx)/2))
    plt.yticks(np.arange(miny, maxy+maxy/2, (maxy-miny)/2))
    if param=='rot':
        plt.xlabel(r'${\bf{\xi}}_{\bf{R} x} (rad)$',fontsize=20, labelpad=8)
        plt.ylabel(r'${\bf{\xi}}_{\bf{R} z} (rad)$',fontsize=20, labelpad=-8)
    if param=='trans':
        plt.xlabel(r'${\bf{\xi}}_{\bf{t} x}(mm)$',fontsize=20, labelpad=8)
        plt.ylabel(r'${\bf{\xi}}_{\bf{t} z}(mm)$',fontsize=20, labelpad=-8)
    ax.set(aspect='equal')
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.8, hspace=None)
    # plt.show(True)
    
    return True


def FCParkSolution(alpha,beta,ta,tb):
    # FCPark solution
    # RotX
    M = np.zeros(shape=(3,3))
    for j in range(len(alpha)):
        M = M + np.asmatrix(beta[j].reshape((3,1)))*np.asmatrix(alpha[j].reshape((3,1))).T
        eig_val,eig_vec = np.linalg.eig(M.T * M)
    FCPark_Rx =  np.asarray(eig_vec*np.diag(np.sqrt(1.0/eig_val))*np.linalg.inv(eig_vec)*M.T)
    # Estimate tx
    C = np.eye(3)-SE3.VecToRot(alpha[0])
    for i in range(1,len(alpha)):
        C = np.vstack((C,np.eye(3)-SE3.VecToRot(alpha[i])))
    g = ta[0] - np.dot(FCPark_Rx,tb[0])
    for i in range(1,len(alpha)):
        g = np.vstack((g, ta[i] - np.dot(FCPark_Rx,tb[i])))
    g = g.reshape(3*len(alpha),1)
    FCPark_tx = np.dot(np.linalg.pinv(C),g).reshape(3)
    return FCPark_Rx, FCPark_tx

def IterativeSolutionTrans(beta, alpha, ta, tb, Rx, sigmaRa, sigmaRb, sigmata, sigmatb, sigmaRx,sigmaRbeta,txinit=np.zeros((3,1)), max_iter =10):
    tx = txinit.reshape(3, 1)
    i = 0
    delta_tx = np.ones((3,1))
    delta_xiRak = np.ones((3,1))
    Ra = []
    Rb = []
    # compute covariance
    inv_sigmaX = []
    for k in range(len(alpha)):
        Ra.append(SE3.VecToRot(alpha[k]))
        Rb.append(SE3.VecToRot(beta[k]))
        sigmaXk = np.zeros((6,6))
        sigmaXk[:3,:3] = sigmaRa[k]
        Rxtbk = np.dot(Rx,tb[k])
        sigmaXk[3:6,3:6] = sigmata + np.dot(np.dot(Rx,sigmatb),np.transpose(Rx)) + np.dot(np.dot(SE3.Hat(Rxtbk),sigmaRx),np.transpose(SE3.Hat(Rxtbk))) #sigmaqk
        transposeJacalpha = np.transpose(SE3.VecToJac(alpha[k]))
        sigmaRxRa = np.dot(np.dot(sigmaRx,np.transpose(SE3.Hat(np.dot(Rx,beta[k])))),transposeJacalpha )-np.dot(sigmaRbeta[k],np.dot(np.transpose(Rx),transposeJacalpha))
        sigmaqRa= -np.dot(SE3.Hat(np.dot(Rx,ta[k])),sigmaRxRa)
        sigmaXk[3:,:3] = sigmaqRa
        sigmaXk[:3,3:] = np.transpose(sigmaqRa)
        inv_sigmaX.append(np.linalg.inv(sigmaXk))
    # import IPython; IPython.embed()
    Rahat = copy.copy(Ra)
    # main loop
    while np.linalg.norm(delta_tx) > 1e-5 and i < max_iter:
        Ak = np.zeros((6,3))
        Bk = np.zeros((6,3))
        Bk[:3,:] = np.eye(3)
        U = np.zeros((3,3))
        eA = np.zeros((3,1))
        sum1 = np.zeros((3,3))
        sum2 = np.zeros((3,1))
        V = []
        W = []
        eB = []
        qhat = []
        [qhat.append(np.dot(Ra_k-np.eye(3),tx)) for Ra_k in Rahat]
        for k in range(len(alpha)):
            Ak[3:6,:] = Ra[k]-np.eye(3)
            # print "Rahat[k]", Rahat[k]
            # print "tx", tx
            Bk[3:6,:] = - SE3.Hat(np.dot(Rahat[k], tx))
            U+= np.dot(np.dot(np.transpose(Ak),inv_sigmaX[k]),Ak)
            Vk = np.dot(np.dot(np.transpose(Bk),inv_sigmaX[k]),Bk)
            Wk = np.dot(np.dot(np.transpose(Ak),inv_sigmaX[k]),Bk)
            V.append(Vk)
            W.append(Wk)
            Yk = np.dot(Wk,np.linalg.inv(Vk))
            Xk = np.zeros((6,1))
            Xk[:3] = np.zeros((3,1))
            Xk[3:6] = (np.dot(Rx,tb[k])-ta[k]).reshape(3,1)
            Xkhat = np.zeros((6,1))
            if math.isnan((SE3.RotToVec(np.dot(Rahat[k], np.linalg.inv(np.array(Ra[k], dtype='float'))))).reshape((3,1))[0]):
                Xkhat[:3] = np.zeros((3,1))
            else:
                Xkhat[:3] = (SE3.RotToVec(np.dot(Rahat[k], np.linalg.inv(np.array(Ra[k], dtype='float'))))).reshape((3,1))
            Xkhat[3:6] = qhat[k].reshape(3,1)
            ek = Xk - Xkhat
            eA += np.dot(np.dot(np.transpose(Ak),inv_sigmaX[k]),ek)
            eBk = np.dot(np.dot(np.transpose(Bk),inv_sigmaX[k]),ek)
            eB.append(eBk)
            sum1 += np.dot(Yk,np.transpose(Wk))
            sum2 += np.dot(Yk,eBk)
        delta_tx = np.dot(np.linalg.inv(U - sum1),eA - sum2)
        tx = tx + delta_tx

        for k in range(len(alpha)):
            delta_xiRak = np.dot(np.linalg.inv(V[k]),(eB[k] - np.dot(np.transpose(W[k]),delta_tx)))
            Rahat[k] = np.dot(SE3.VecToRot(delta_xiRak),Rahat[k])
        i = i+1
    # sigmatx = np.linalg.inv(U - sum1)
    if i < max_iter:
      return tx, np.linalg.inv(U - sum1), i
    else:
        return tx, np.linalg.inv(U - sum1), False

def IterativeSolutionRot(beta,alpha,sigmaRa,sigmaRb,Rxinit=np.eye(3),max_iter = 10):
    # Errors in A and B
    Rhat = Rxinit
    i = 0
    betahat = copy.copy(beta)
    deltaA = np.ones(3)
    deltaBk = np.ones(3)
    # TODO: changing the names of the variables
    inv_sigmaX= []
    for k in range(len(alpha)):
        sigmaXk = np.zeros((6,6))
        sigmaXk[3:6,3:6] = np.dot(np.dot(SE3.VecToJacInv(alpha[k]),sigmaRa),np.transpose(SE3.VecToJacInv(alpha[k])))
        sigmaXk[:3,:3] = np.dot(np.dot(SE3.VecToJacInv(beta[k]),sigmaRb),np.transpose(SE3.VecToJacInv(beta[k])))
        inv_sigmaXk = np.linalg.inv(sigmaXk)
        inv_sigmaX.append(inv_sigmaXk)
    while np.linalg.norm(deltaA)>1e-10 and i < max_iter:
        Ak = np.zeros((6,3))
        Bk = np.zeros((6,3))
        Bk[:3,:] = np.eye(3)
        U = np.zeros((3,3))
        eA = np.zeros((3,1))
        sum1 = np.zeros((3,3))
        sum2 = np.zeros((3,1))
        V = []
        W = []
        eB = []
        alphahat = []
        Y  = []
        [alphahat.append(np.dot(Rhat,betak)) for betak in betahat]
        for k in range(len(alpha)):
            Ak[3:6,:] = -SE3.Hat(np.dot(Rhat,betahat[k]))
            Bk[3:6,:] = Rhat
            U += np.dot(np.dot(np.transpose(Ak),inv_sigmaX[k]),Ak)
            Vk = np.dot(np.dot(np.transpose(Bk),inv_sigmaX[k]),Bk)
            Wk = np.dot(np.dot(np.transpose(Ak),inv_sigmaX[k]),Bk)
            V.append(Vk)
            W.append(Wk)
            Yk = np.dot(Wk,np.linalg.inv(Vk))
            Y.append(Yk)
            Xk = np.zeros((6,1))
            Xk[:3] = beta[k].reshape(3,1)
            Xk[3:6] = alpha[k].reshape(3,1)
            Xkhat = np.zeros((6,1))
            Xkhat[:3] = betahat[k].reshape(3,1)
            Xkhat[3:6] = alphahat[k].reshape(3,1)
            ek = Xk - Xkhat
            eA += np.dot(np.dot(np.transpose(Ak),inv_sigmaX[k]),ek)
            eBk = np.dot(np.dot(np.transpose(Bk),inv_sigmaX[k]),ek)
            eB.append(eBk)
            sum1 += np.dot(Yk,np.transpose(Wk))
            sum2 += np.dot(Yk,eBk)
        deltaA = np.dot(np.linalg.inv(U - sum1),eA - sum2)
        Rhat = np.dot(SE3.VecToRot(deltaA), Rhat)
        for k in range(len(alpha)):
            deltaBk = np.dot(np.linalg.inv(V[k]),(eB[k] - np.dot(np.transpose(W[k]),deltaA)))
            betahat[k] = betahat[k] + deltaBk.reshape(3)
        i +=1
    if i < max_iter:
        # import IPython; IPython.embed()
        alphahat = []
        [alphahat.append(np.dot(Rhat,betak)) for betak in betahat]
        sigmaRx = np.linalg.inv(U - sum1)
        sigmaRbeta = [-np.dot(sigmaRx,Yn) for Yn in Y]
        sigmabeta = [np.dot(np.dot(np.transpose(Y[n]),sigmaRx),Y[n])+np.linalg.inv(V[n]) for n in range(len(alpha))]
        sigmanewRa = []
        sigmaRRa = []
        for n in range(len(alphahat)):
            Jacalpha = SE3.VecToJac(alphahat[n])
            JachatRbeta = -np.dot(Jacalpha,SE3.Hat(np.dot(Rhat,betahat[n])))
            JacR = np.dot(Jacalpha,Rhat)
            sigmanewRa.append(np.dot(JachatRbeta,np.dot(sigmaRx,np.transpose(JachatRbeta)))+np.dot(JacR,np.dot(sigmabeta[n],np.transpose(JacR)))+np.dot(JacR,np.dot(np.transpose(sigmaRbeta[n]),np.transpose(JachatRbeta)))+np.dot(JachatRbeta,np.dot(sigmaRbeta[n],np.transpose(JacR))))
            sigmaRRa.append(np.dot(sigmaRbeta[n],np.transpose(JacR))+np.dot(sigmaRx,np.transpose(JachatRbeta)))
        return Rhat, sigmaRx,i,betahat,alphahat,sigmaRbeta,sigmabeta,sigmanewRa,sigmaRRa
    else:
        return Rhat, np.linalg.inv(U - sum1), False,beta,alpha,None,None,None,None
