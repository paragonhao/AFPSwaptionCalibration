#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 15:31:08 2019

@author: paragonhao
"""

class G2PlusPlusParams:
    
    alpha = None 
    beta = None 
    eta = None 
    sigma = None 
    rho = None 
    
    def __init__(self, alpha, beta, eta, sigma, rho):
        self.alpha = alpha
        self.beta = beta 
        self.eta = eta 
        self.sigma = sigma 
        self.rho = rho 
        
    def printParams(self):
        print("The Parameters values are: ")
        print("alpha: {}".format(self.alpha))
        print("beta: {}".format(self.beta))
        print("eta: {}".format(self.eta))
        print("sigma: {}".format(self.sigma))
        print("rho: {}".format(self.rho))