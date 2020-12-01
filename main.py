# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 16:57:59 2020

@author: Andreas
"""

import numpy as np
import scipy.linalg as la
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
# from plotting import ellipse
#from utils import rotmat2d

from tqdm import tqdm
from plot_setup import setup_plot
from plotting import *
from plotting import *
setup_plot()

# %% Ship parameters
m = 17.0677e6                      # mass (kg)
Iz = 2.1732e10                     # yaw moment of inertia about CO (kg m^3)
xg = -3.7                          # CG x-ccordinate (m)
L = 161                            # length (m)
B = 21.8                           # beam (m)
T = 8.9                            # draft (m)
Dia = 3.3                          # propeller diameter (m)
rho = 1025                         # density of water (kg/m^3)
visc = 1e-6                        # kinematic viscousity at 20 degrees (m/s^2)
eps = 0.001                        # a small number added to ensure that the denominator of Cf is well defined at u=0
k = 0.1                            # form factor giving a viscous correction
t_thr = 0.05                       # thrust deduction number

# Get KT and KQ from Wagening B-series propeller table
J_a = 0                            # Bollard Pull

AEAO = 0.65
PD = 1.5                           # Pitch diameter ratio
z = 4                              # Number of blades

# Rudder limitations
delta_max = 40 * np.pi/180
Ddelta_max = 5 * np.pi/180

# Added mass matrix about CO
Xudot = -8.9830e5
Yvdot = -5.1996e6
Yrdot =  9.3677e5
Nvdot =  Yrdot
Nrdot = -2.4283e10

MA = -np.array([
    [Xudot, 0, 0],
    [0, Yvdot, Yrdot],
    [0, Nvdot, Nrdot]
    ])

# Rigid body mass matrix
MRB = np.array([
    [m, 0, 0],
    [0, m, m*xg],
    [0, m*xg, Iz]
    ])

#Total inertia mass matrix inverted:
Minv = la.inv(MRB + MA)

# %% Environmental disturbance models
# The model for current and wind is implemented here

# %% Ocean current in NED
Vc = 0                          # Current speed (m/s)
betaVc = np.deg2rad(45)         #Current direction (rad)
nu_c = np.transpose(
    np.array([
    [-np.cos(betaVc)*Vc, -np.sin(betaVc)*Vc, 0]
    ])
    )

# %% Wind expressed in NED
Wv = 0                          # Wind speed (m/s)
betaWv = deg2rad(135)           # Wind direction (rad)
rho_a = 1.247                   # Air density at 10 deg celsius
cy = 0.95                       # wind coefficient in sway
cn = 0.15                       # wind coefficient in yaw
A_lw = 10 * L                   # projected lateral area of ownship