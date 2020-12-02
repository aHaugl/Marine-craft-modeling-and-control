# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 19:04:43 2020

@author: Andreas
"""

from typing import Tuple
import numpy as np
from scipy.linalg import block_diag
import scipy.linalg as la
from utilities import rotmat2d
from JCBB import JCBB
import utils
from utils import wrapToPi
from numpy import matlib as ml


class Controller:
    def __init__(
    self,
    Q,
    R,
    do_asso=False,
    alphas=np.array([0.001, 0.0001]),
    sensor_offset=np.zeros(2),
):

    self.Q = Q
    self.R = R
    self.do_asso = do_asso
    self.alphas = alphas
    self.sensor_offset = sensor_offset
    
    def 
