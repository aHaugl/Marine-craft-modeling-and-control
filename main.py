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
from utilities import *
from plot_setup import setup_plot
from plotting import *
from plotting import *
setup_plot()

# %% User inputs
h = 0.1                            # Sampling time [s]
Ns = 10000                         # Number of samples
loadmat("WP.mat")                  # Waypoints

psi_ref = 0 * np.pi / 180          # desired yaw angle (rad)
U_d = 7                            # desired cruise speed (m/s)

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

# initial states
eta = np.array([[0,    0,  0]]).T
nu  = np.array([[0.1,  0,  0.1]]).T
delta = 0;
n = 0;

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

# Linear damping matrix (only valid for zero speed)
T1 = 20
T2 = 20
T6 = 10

Xu = -(m - Xudot) / T1
Yv = -(m - Yvdot) / T2
Nr = -(Iz - Nrdot) / T6
D = np.diag([-Xu, -Yv, -Nr])

# Get KT and KQ from Wagening B-series propeller table
J_a = 0                            # Bollard Pull
AEAO = 0.65
PD = 1.5                           # Pitch diameter ratio
z = 4                              # Number of blades

[KT, KQ ] = Wageningen(J_a, PD, AEAO, z)


# Rudder limitations
delta_max = 40 * np.pi/180
Ddelta_max = 5 * np.pi/180

# Rudder coefficients (Section 9.5 Fossen2021)
b = 2
AR = 8
CB = 0.

Lambda = b**2 / AR
tR = 0.45 - 0.28 * CB
CN = 6.13 * Lambda / (Lambda + 2.25)
aH = 0.75
xH = -0.4 * L
xR = -0.5 * L

X_delta2 = 0.5 * (1 - tR) * rho * AR * CN
Y_delta = 0.25 * (1 + aH) * rho * AR * CN;
N_delta = 0.25 * (xR + aH * xH) * rho * AR * CN

# input matrix
t_thr = 0.05;           # thrust deduction number
X_delta2 = 0;           # rudder coefficients (Section 9.5)
Y_delta = 0;      
N_delta = 1;
Bu = np.array([
    [(1 - t_thr),    X_delta2],
    [0,              Y_delta],
    [0,              N_delta]
    ])

# %% Environmental disturbance models
# The model for current and wind is implemented here

# %% Ocean current in NED
Vc = 0                          # Current speed (m/s)
betaVc = np.deg2rad(45)         #Current direction (rad)
nu_c =  np.array([
    [-np.cos(betaVc)*Vc, -np.sin(betaVc)*Vc, 0]
    ]).T

# %% Wind expressed in NED
Wv = 0                          # Wind speed (m/s)
betaWv = np.deg2rad(135)           # Wind direction (rad)
rho_a = 1.247                   # Air density at 10 deg celsius
cy = 0.95                       # wind coefficient in sway
cn = 0.15                       # wind coefficient in yaw
A_lw = 10 * L                   # projected lateral area of ownship

# %% A wave model can be expressed here, but were not required for the project

# %% MAIN LOOP


simdata = np.zeros((Ns, 18))                # table of simulation data

controller_heading = np.zeros((Ns+1, 5))            #table of controller data

for i in range (Ns):

    t = (i-1) * h;                        # time (s)
   
    # %% state-dependent time-varying matrices
    CRB = m * nu(3) * np.array([
                            [0, -1, -xg ],
                            [1,  0,  0   ],
                            [xg, 0,  0  ]
                            ])
    
    R = Rzyx(0,0,eta(3));
    
    # %% reference models
    psi_d = psi_ref
    r_d = 0
    u_d = u_ref
   
    # %% thrust 
    thr = rho * Dia^4 * KT * abs(n) * n        # thrust command (N)
        
    # %% control law
    delta_c = 0.1                              # rudder angle command (rad)
    n_c = 10                                   # propeller speed (rps)
    
    # %% ship dynamics
    u = np.array([ thr, delta ]).T
    tau = B * u
    nu_dot = Minv * (tau - CRB * nu)
    eta_dot = R * nu    
    
    # Rudder saturation and dynamics (Sections 9.5.2)
    if abs(delta_c) >= delta_max:
        delta_c = sign(delta_c)*delta_max;
    end
    
    delta_dot = delta_c - delta;
    if abs(delta_dot) >= Ddelta_max:
        delta_dot = sign(delta_dot)*Ddelta_max;
    end    
    
    # propeller dynamics
    n_dot = (1/10) * (n_c - n);
    
    # %% store simulation data in a table (for testing)
    simdata(i,:) = np.concatenate(t, n_c, delta_c, n, delta, eta.T, nu.T,
                              u_d, psi_d, r_d, psi_meas, x_hat.T, r_meas)    
     
    # %% Euler integration
    eta = euler2(eta_dot,eta,h); 
    nu  = euler2(nu_dot,nu,h);
    delta  = euler2(delta_dot,delta,h);   
    n  = euler2(n_dot,n,h);    
    
end
