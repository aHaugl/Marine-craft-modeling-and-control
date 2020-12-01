# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 17:30:14 2020

@author: Andreas
"""

import numpy as np

def wrapToPi(angle):
    return (angle + np.pi) % (2*np.pi) - np.pi

def rotmat2d(angle):
    return np.array([[np.cos(angle), -np.sin(angle)],
                     [np.sin(angle),  np.cos(angle)]])
