#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 15:51:31 2019

@author: paragonhao
"""
import math 

# G2++ zero coupon bond analytical formula
class SWPTNG2PPAF:
    
    termStructure = None
    
    def __init__(self, termStructure):
        self.termStructure = termStructure
        termStructure['Term[y]'] = termStructure['Term[y]'].apply(lambda x: round(float(x),1))
        
    
    
    def getTermStructure(self, T):
        currTS = self.termStructure[self.termStructure['Term[y]']==T]
        
        return float(currTS['DF.mid'].mean())
    
    
    
    # V(t,T) in the book page 145
    # g2params: the five parameters 
    # t: starting time of the ZCB 
    # T: ending time of the ZCB
    def V(self, g2params, t):
        # five params in an array 
        alpha = g2params[0]
        beta = g2params[1]
        sigma = g2params[2]
        eta = g2params[3]
        rho = g2params[4]
        
        expat = math.exp(-alpha * t)
        expbt = math.exp(-beta * t)
        cx = sigma/alpha
        cy = eta/beta
        temp1 = cx * cx * (t + (2.0 * expat - 0.5 * expat * expat - 1.5)/alpha)
        temp2 = cy * cy * (t + (2.0 * expbt - 0.5 * expbt * expbt - 1.5)/beta)
        temp3 = 2.0 * rho * cx * cy * (t + (expat -1.0)/alpha + (expbt - 1.0)/beta - 
                                       (expat * expbt - 1.0)/(alpha + beta))
        
        return temp1 + temp2 + temp3
        
    
    
    # A(t, T) function in the book page 148 
    # g2params: the five parameters 
    # Note: t_i > T
    def A(self, g2params, t_i, T):
        return (self.getTermStructure(T)/self.getTermStructure(t_i)) * math.exp(0.5 * (self.V(g2params, T - t_i) - 
                                      self.V(g2params, T) + self.V(g2params, t_i)))


    
    def B(self, x, t):
        return (1.0 - math.exp(-x*t))/x
        
    
    
    # mu_x function in book page 154
    # mu 1 
    def Mu_x(self, g2params, T):
        # five params in an array 
        alpha = g2params[0]
        beta = g2params[1]
        sigma = g2params[2]
        eta = g2params[3]
        rho = g2params[4]
        
        temp = sigma * sigma/(alpha * alpha)
        
        mux = -(temp + rho * sigma * eta/(alpha * beta)) * (1.0 - math.exp(-alpha * T)) \
        - 0.5 * temp * (1.0 - math.exp(-2.0 * alpha * T )) - rho * sigma * eta/(beta * (alpha + beta)) \
        * (1 - math.exp(-(beta + alpha) * T))
        
        return mux 



    # mu_y function in book page 154
    # mu 2
    def Mu_y(self, g2params, T):
        # five params in an array 
        alpha = g2params[0]
        beta = g2params[1]
        sigma = g2params[2]
        eta = g2params[3]
        rho = g2params[4]
        
        temp = eta * eta /(beta * beta)
        
        muy = -(temp + rho * sigma * eta/(alpha * beta)) * (1.0 - math.exp(-beta * T)) \
        - 0.5 * temp * (1.0 - math.exp(-2.0 * beta * T)) \
        - rho * sigma * eta /(alpha * (alpha + beta)) * (1.0 - math.exp(-(beta + alpha) * T))
        
        return muy
    



    def sigma_x(self, sigma, alpha, T):
        sigma_x = sigma * math.sqrt(0.5 * (1.0 - math.exp(-2.0 * alpha * T))/alpha)
        return sigma_x
    
    
    
    def sigma_y(self, eta, beta, T):
        sigma_y = eta * math.sqrt(0.5 * (1.0 - math.exp(-2.0 * beta * T))/beta)
        return sigma_y
    
    

    def rho_xy(self, g2params, sigma_x, sigma_y, T):
        alpha = g2params[0]
        beta = g2params[1]
        sigma = g2params[2]
        eta = g2params[3]
        rho = g2params[4]
        
        return rho * eta * sigma * (1.0 - math.exp(-(alpha + beta) * T)) /((alpha + beta) * sigma_x * sigma_y)
         
    
    