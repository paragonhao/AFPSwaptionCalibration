#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 12:19:03 2019

@author: paragonhao
"""
import math
from scipy.stats import norm


class BlackModel76:
    
    tenor = None  #tenor of swap in years
    fRate = None  #forward rate of underlying swap 
    xRate = None  #strike rate of swaption
    rfRate = None #risk-free interest rate
    T = None      #time to expiration in years
    sigma = None  #volatility of the forward-starting swap rate
    m = None      #compoundings per year in swap rate
    d1 = None
    d2 = None
    df = None
    
    
    def __init__(self, tenor, fRate, xRate, rfRate, T, sigma, m):
        self.tenor = tenor 
        self.fRate = fRate 
        self.xRate = xRate
        self.rfRate = rfRate
        self.T = T 
        self.sigma = sigma 
        self.m = m 
        self.d1 = (math.log(fRate/xRate) + 0.5 * (sigma ** 2) * T)/(sigma * math.sqrt(T))
        self.d2 = self.d1 - sigma * math.sqrt(T)
        self.df = math.exp(-1 * self.rfRate * self.T)
    
    
    
    def getPayerSwaptionPrice(self):
        temp = 1/ (1 + (self.fRate/self.m)) ** (self.tenor * self.m)
        nd1 = norm.cdf(self.d1)
        nd2 = norm.cdf(self.d2)
        return ((1 - temp)/self.fRate) * self.df * (self.fRate * norm.cdf(self.d1) -  self.xRate * norm.cdf(self.d2))
        
    
    
    def getReceiverSwaptionPrice(self):
        temp = (1 + (self.fRate/self.m)) ** (self.tenor * self.m)
        nd1_neg = norm.cdf(-1 * self.d1)
        nd2_neg = norm.cdf(-1 * self.d2)
        return ((1 - temp)/self.fRate) * self.df * (self.xRate * nd2_neg - self.fRate * nd1_neg)
    