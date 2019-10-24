#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 20:32:45 2019

@author: paragonhao
"""
import math
from scipy.stats import norm

class NormalModel:
    
    fRate = None  #forward rate of underlying swap 
    strike = None
    impliedVol = None
    maturity = None
    d1 = None
    discountFactors = None
    
    
    def __init__(self, fRate, strike, impliedVol, maturity, discountFactors):
        self.fRate = fRate
        self.strike = strike
        self.impliedVol = impliedVol
        self.maturity = maturity
        self.discountFactors = discountFactors
        self.d1 = (fRate - strike)/(impliedVol * math.sqrt(maturity))
        
    def getPayerSwaptionPrice(self):
        temp = (self.fRate - self.strike) * norm.cdf(self.d1) + self.impliedVol * math.sqrt(self.maturity) * norm.pdf(self.d1)
        return sum(self.discountFactors) * temp
    
    
    def getReceiverSwaptionPrice(self):
        temp = (self.strike - self.fRate) * norm.pdf(-1 * self.d1) + self.impliedVol * math.sqrt(self.maturity) * norm.cdf(self.d1)
        return sum(self.discountFactors) * temp