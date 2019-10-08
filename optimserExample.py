#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 15:04:22 2019

@author: paragonhao
"""

import scipy.optimize as optimize
import math


def f(params):
    # print(params)  # <-- you'll see that params is a NumPy array
    a,b = params # <-- for readability you may wish to assign names to the component variables
    return abs(a + b**2 - 4)

initial_guess = [2, 1]
bds = ((0, None),(0, None))

result = optimize.minimize(f, initial_guess, method='Nelder-Mead')

if result.success:
    fitted_params = result.x
    print(fitted_params)
else:
    raise ValueError(result.message)
    
from scipy.optimize import newton 

def f(x):
    return x**2 - 4


print(newton(f, 5))