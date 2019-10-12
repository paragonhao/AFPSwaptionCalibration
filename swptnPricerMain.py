#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 16:59:00 2019

@author: paragonhao
"""


from scipy.stats import norm
# norm.cdf(0)
import numpy as np
from swptnG2PPAF import SWPTNG2PPAF
from utils import Utils
import math
import scipy.optimize as optimize

import datetime
import pandas as pd

# parametmers of the swaption to calibrate  
# tenor of the IRS
tenor = 1.0

# the maturity of the option on the swap
maturity = 2.0

# total notional value for the IRS
notional = 1.0

# fixed rate
fixedRate = 0.0166

# determine if it is a payer or receiver swaption
isPayer = 1

# frequency of the payment 
payFreq = 0.5

#load DFfactors
termStructure = Utils.readDF('termStructure.csv')

G2 = SWPTNG2PPAF(termStructure)

def swaptionPricingFunction(g2params, tenor, maturity, notional, fixedRate, isPayer):
    print("running the swaptionPricingFunction")
    # five params in an array 
    alpha = g2params[0]
    beta = g2params[1]
    sigma = g2params[2]
    eta = g2params[3]
    rho = g2params[4]
    print("paramss+++++")
    print(g2params)
    print("paramss+++++")
    P_0_T = G2.getTermStructure(maturity)
    
    paymentTimes = tenor/payFreq
    tau = np.repeat(payFreq, paymentTimes, axis=0)
    
    c_i = tau * fixedRate
    c_i[-1] = c_i[-1] + 1
    t_i = np.arange(maturity+payFreq, tenor+maturity + payFreq, payFreq)
#    print("paymentTimes: {}".format(len(t_i)))
#    print("cash flow: {}".format(c_i))
          
    sigma_x = G2.sigma_x(sigma, alpha, maturity)
    
    sigma_y = G2.sigma_y(eta, beta, maturity)
    
    mu_x = G2.Mu_x(g2params, maturity)
    
    mu_y = G2.Mu_y(g2params, maturity)
    
    rho_xy = G2.rho_xy(g2params, sigma_x, sigma_y, maturity)
    
    txy = math.sqrt(1 - rho_xy * rho_xy)
    
#    print("sigma_x: {}".format(sigma_x))
#    print("sigma_y: {}".format(sigma_y))
#    print("mu_x: {}".format(mu_x))
#    print("mu_y: {}".format(mu_y))
#    print("rho_xy: {}".format(rho_xy))
#    print("txy: {}".format(rho_xy))
    
    A_array = []
    Ba_array = []
    Bb_array = []
    lambda_array = np.repeat(0, len(t_i), axis=0)
    
    # define the boundary 
    interval = 0.01
    lower = mu_x - 10 * sigma_x
    upper = mu_x + 10 * sigma_x
    xList = np.arange(lower, upper, interval)
    
#    print("Upper: {}".format(upper))
#    print("Lower: {}".format(lower))
    
    for i in range(len(t_i)):
#        print("range: {}".format(i))
        A_array.append(G2.A(g2params, maturity, t_i[i]))
#        print("A: {}".format(G2.A(g2params, maturity, t_i[i])))
        Ba_array.append(G2.B(alpha, t_i[i] - maturity))
#        print("Ba: {}".format(G2.B(alpha, t_i[i] - maturity)))
        Bb_array.append(G2.B(beta, t_i[i] - maturity))
    
#    print(A_array)
#    print(Ba_array)
#    print(Bb_array)
    
    integral_result = 0 
    
    for x in xList:
        lambda_array = []
        for i in range(len(t_i)):
            lambda_array.append(c_i[i] * A_array[i] * math.exp(-Ba_array[i] * x))

        def y_bar_solver(y):
            sum_all = 0
            
            for i in range(len(t_i)):
                sum_all += lambda_array[i] * math.exp(-Bb_array[i] * y)
            return sum_all - 1
        
        y_bar = optimize.fsolve(y_bar_solver, x0=1)
        h1 = (y_bar - mu_y)/(sigma_y * txy) - rho_xy * (x - mu_x)/(sigma_x * txy)
        cdf_val_1 = norm.cdf(-isPayer * h1)
        
        try:
            for i in range(len(t_i)):
#                print("Bb_array[i]: {}".format(Bb_array[i]))
#                print("lambda_array[i]: {}".format(lambda_array[i]))
                h2 = h1 + Bb_array[i] * sigma_y * txy
                kappa = - Bb_array[i] * (mu_y - 0.5 * txy * txy * sigma_y * sigma_y * Bb_array[i] \
                                  + rho_xy * sigma_y * (x - mu_x)/sigma_x)
#                print("kappa: {}".format(kappa))
                cdf_val_1 -= lambda_array[i] * math.exp(kappa) * norm.cdf(-h2 * isPayer)[0]
        except OverflowError:
            print("Overflow err: ")
            
        temp_pdf = (x - mu_x) / sigma_x
        
        integral_result += math.exp(-0.5 * temp_pdf * temp_pdf) * cdf_val_1/(sigma_x * math.sqrt(2.0 * math.pi))
    
    swaptionPrice = notional * isPayer * P_0_T * integral_result
    
    
    # getting implied volatilty
    c = swaptionPrice /(isPayer * notional * (G2.getTermStructure(maturity) - G2.getTermStructure(maturity + tenor)))
    g2_IV = (1/(isPayer * math.sqrt(maturity))) * (2.506297 * c - 0.686461 * c * c)/(1 - 0.277069 * c - 0.237552 * c * c)
    return g2_IV - 0.709
    
g2params = [2.8187,0.035,0.0579,0.0091,-0.999]
swaptionPricingFunction(g2params, tenor, maturity, notional, fixedRate, isPayer)
bnds = ((0.001, 5),(0.001, 0.05),(0.001, 0.05),(0.001, 0.05),(-0.999, 0.999))

result = optimize.minimize(swaptionPricingFunction, g2params, args=(tenor, maturity, notional, fixedRate, isPayer), bounds=bnds, method='L-BFGS-B')

if result.success:
    fitted_params = result.x
    print("current fitted value")
    print(fitted_params)
else:
    raise ValueError(result.message)






