#!/usr/bin/env python3
# Christopher Arran
# Created 4/9/19
# A very quick interferometry analysis toolbox

import scipy.constants as constant
import numpy

def phaseToDensity(phase,L,lambda0):    
    omega0 = 2*constant.pi*constant.c/lambda0
    omegap_sq = (2*omega0*constant.c)*phase/L
    ne = (constant.m_e*constant.epsilon_0/constant.e**2) * omegap_sq
    return ne