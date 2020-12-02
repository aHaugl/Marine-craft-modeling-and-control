# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 17:30:14 2020

@author: Andreas
"""

import numpy as np
from scipy.io import loadmat

def wrapToPi(angle):
    return (angle + np.pi) % (2*np.pi) - np.pi

def rotmat2d(angle):
    return np.array([[np.cos(angle), -np.sin(angle)],
                     [np.sin(angle),  np.cos(angle)]])


# %% MSS toolbox functions from Matlab 
def Rzyx(phi, theta, psi):
    """
    # R = Rzyx(phi,theta,psi) computes the Euler angle
    # rotation matrix R in SO(3) using the zyx convention
    
    # Author:       Thor I. Fossen
    # Date:         14th June 2001
    # Revisions:    Translation to Python
    # Author:       Andreas Haugland
    # Date:         02.12.2020
    """
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    cth  = np.cos(theta)
    sth  = np.sin(theta)
    cpsi = np.cos(psi)
    spsi = np.sin(psi)
    
    return np.array([
        [cpsi*cth,      -spsi*cphi+cpsi*sth*sphi,    spsi*sphi+cpsi*cphi*sth],
        [spsi*cth,      cpsi*cphi+sphi*sth*spsi,     -cpsi*sphi+sth*spsi*cphi],
        [-sth,          cth*sphi,                    cth*cphi]
        ])

def Wageningen(J_a, PD, AEAO, z):
    """ [KT, KQ] = wageningen(Ja,PD,AEAO,z) computes the thrust and torque
     coefficients of the Wageningen B-series propellers for:
    
     Ja = Va/(n*D)  advance number
     PD             pitch/diameter ratio (typically 0.5-2.5)
     AEAO           blade area ratio (ratio of expanded blade area to propeller disk area)
     z              number of propeller blades
    
     see ExWageningen.m for use of [K_T, K_Q] = Wageningen(Ja,PD,AEAO,z)
    
     The B-series propellers were designed and tested at the Netherlands Ship
     Model Basin in Wageningen. The open_water characteristics of 120
     propeller models of the B_series were tested at N.S.M.B. and analyzed
     with multiple polynomial regression analysis. The derived polynomials
     express the thrust and torque coefficients in terms of the number of
     blades, the blade area ratio, the pitch_diameter ratio and the advance
     coefficient.
    
     Reference: Barnitsas, M.M., Ray, D. and Kinley, P. (1981).
     KT, KQ and Efficiency Curves for the Wageningen B-Series Propellers
     http://deepblue.lib.umich.edu/handle/2027.42/3557
    
     Author:    Thor I. Fossen
     Date:      7 October 2018
     Revisions: 18 June 2019 - minor updates
     
     Translation to Python:
     Author:    Andreas Haugland
     Date:      02.12.2020
    """
    
    WageningenMat = loadmat("WageningData.mat")
    WagCThrust_stuv = WageningenMat["WagCThrust_stuv"]
    WagThrust_s = WageningenMat["WagThrust_s"]
    WagThrust_t = WageningenMat["WagThrust_t"]
    WagThrust_u = WageningenMat["WagThrust_u"]
    WagThrust_v = WageningenMat["WagThrust_v"]
    
    WagCTorque_stuv = WageningenMat["WagCTorque_stuv"]
    WagTorque_s = WageningenMat["WagTorque_s"]
    WagTorque_t = WageningenMat["WagTorque_t"]
    WagTorque_u = WageningenMat["WagTorque_u"]
    WagTorque_v = WageningenMat["WagTorque_v"]
    
    KT = sum(WagCThrust_stuv * (J_a ** WagThrust_s) * (PD ** WagThrust_t) 
             * (AEAO ** WagThrust_u) * (z ** WagThrust_v)
              )
    KQ = sum(WagCTorque_stuv * (J_a ** WagTorque_s) * (PD ** WagTorque_t)
             * (AEAO ** WagTorque_u) * (z ** WagTorque_v)
             )

    return np.array([KT, KQ])