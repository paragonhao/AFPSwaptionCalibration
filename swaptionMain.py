#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 15:27:24 2019

@author: paragonhao
"""

from scipy.stats import norm
# norm.cdf(0)
import numpy as np
from g2PlusPlusParams import G2PlusPlusParams as G2Parameters
from g2PlusPlusAnalyticalFormula import G2PlusPlusAnalyticalFormula as GA
# example 

g2params = G2Parameters(1,1.5,2,0.16,0.5)
tenor = 1
notional = 100 
X = 0.02
isPayer = 1
maturity = 2
payFreq = 0.5


def swaptionPrice(g2params, T, maturity, notional, X, isPayer):
    
    interval = 0.001
    
    # get the variable c 
    paymentTimes = maturity/payFreq
    yearFractions = np.repeat(payFreq, paymentTimes, axis=0)
    c_i = yearFractions * X
    c_i[-1] = c_i[-1] + 1
    t_i = np.arange(T+payFreq, T+maturity + payFreq, payFreq)
    
    # first part of the equation outside the integral;
    discountedNotional = notional * isPayer * GA.getMarket_DF(T)
    
    sigma_x = GA.sigma_x(g2params, T)
    sigma_y= GA.sigma_y(g2params, T)
    
    mu_x = -1 * GA.Moment_T_x(g2params, 0, T, T)
    mu_y = -1 * GA.Moment_T_y(g2params, 0, T, T)
    
    upper_bound = mu_x + 10 * sigma_x
    lower_bound = mu_x - 10 * sigma_x
    
    rho_xy = GA.rho_xy(g2params, T, sigma_x, sigma_y)
    
    # numerically solve the integration 
    xList = np.arange(lower_bound, upper_bound, interval)
    
    price = 0 
    
    for x in xList:
        
        # first part in the integral
        temp1 = math.exp(-0.5 * (((x - mu_x)/sigma_x) ** 2)) / (sigma_x * math.sqrt( 2 * math.pi))
        
        sum_Func_B = 0 
        sum_Func_Ci_A = 0 
        sum_Func_B_x = 0 
        
        for i in range(len(t_i)):
            sum_Func_B_x += GA.ZCB_Func_B(g2params.alpha, T, t_i[i]) * x
            sum_Func_Ci_A += math.log(c_i[i] * GA.ZCB_Func_A(g2params, T, t_i[i]))
            sum_Func_B += GA.ZCB_Func_B(g2params.alpha, T, t_i[i])
        
        y_bar = (sum_Func_B_x - sum_Func_Ci_A)/sum_Func_B
        h_1_x = ((y_bar - mu_y)/(sigma_y * math.sqrt(1 - rho_xy ** 2))) - ((rho_xy * (x - mu_x))/(sigma_x * math.sqrt(1 - rho_xy ** 2)))
        
        temp2 = norm.cdf(-1 * isPayer * h_1_x)
        
        temp3 = 0
        
        for i in range(len(t_i)):
            lambda_i = c_i[i] * GA.ZCB_Func_A(g2params, T, t_i[i]) * math.exp(-1 * GA.ZCB_Func_B(g2params.alpha, T, t_i[i]) * x)
            kappa_i = -1 * GA.ZCB_Func_B(g2params.beta, T, t_i[i]) * (mu_y - \
                                        0.5 * (1 - rho_xy ** 2) * (sigma_y ** 2) * GA.ZCB_Func_B(g2params.beta, T, t_i[i]) + \
                                        rho_xy * sigma_y * (x - mu_x)/sigma_x)
            
            h_2_x = h_1_x + GA.ZCB_Func_B(g2params.beta, T, t_i[i]) * sigma_y * math.sqrt(1 - rho_xy ** 2)
            
            temp3 += lambda_i * math.exp(kappa_i) * norm.cdf(-1 * isPayer * h_2_x)

        price += temp1 * (temp2 - temp3) * interval
    
    return price 
        

swaptionPrice(g2params, T, maturity, notional, X, isPayer)
        
        
        
        
        
        
    
    